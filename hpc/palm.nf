#!/bin/env nextflow

// Name        : palm.nf
// Description : Fan out independent PALM jobs on an HPC (or locally) with
//               Nextflow + Apptainer/Singularity.
//
// Idea
// ----
// PALM itself is not parallelised. On a cluster you typically run many
// *independent* analyses (different inputs, designs, contrasts, etc.).
// This workflow:
//   1. Reads a shared "template" directory (files staged into every job).
//   2. Reads a task list: one PALM argument line per job.
//   3. Submits one Nextflow process per line.
//   4. Appends a unique -seed and a unique -o to each call.
//
// What belongs in task.list vs template/
// --------------------------------------
//   template/   Shared files every job needs (e.g. mask, surface, design
//               that does not change across tasks). These are staged into
//               the process work directory so relative paths in the task
//               line can refer to them by basename.
//   task.list   One line per independent `palm` invocation. Each line is
//               the CLI arguments *except* -seed and -o (those are added
//               by this workflow). Blank lines and lines starting with '#'
//               are ignored.
//
// Example task.list line:
//   -i mydata.nii -d design.mat -t design.con -m mask.nii -n 5000 -T
//
// Seeds
// -----
// Task k (1-based) uses:  -seed (base_seed + k - 1)
// So with the default --base_seed 1, tasks get seeds 1, 2, 3, ...
// Use a different --base_seed for a later batch if you want that batch's
// RNG streams not to reuse the same integer sequence as an earlier batch.
//
// Outputs
// -------
// Each task writes under out_<k>/ with PALM's -o set to out_<k>/palm, plus
// a small out_<k>/seed.txt recording the seed used. Results are published
// to --dir_out. Optionally pack each task directory into out_<k>.tar.gz.
//
// Typical runs
// ------------
//   nextflow run palm.nf -profile local  --dir_in /path/to/input --dir_out /path/to/out
//   nextflow run palm.nf -profile slurm  --dir_in /path/to/input --base_seed 10001
//
// Input layout expected under --dir_in (unless overridden):
//   dir_in/
//     template/     # shared files
//     task.list     # one PALM arg line per job
//
// Related files
// -------------
//   nextflow.config          executors, container engine, resource labels
//   ../containers/Apptainer.def   builds palm.sif (Octave + PALM)
//
// _______________________________________
// Juan M. Peralta and Anderson M. Winkler
// UTRGV
// June/2026

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PALM -- Permutation Analysis of Linear Models
// Copyright (C) 2015 Anderson M. Winkler
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.
//
// This program is distributed in the hope that it will be useful
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

nextflow.enable.dsl = 2

//============================================================================
// Parameters
// Override on the command line with:  --name value
//============================================================================

// Directory that contains template/ and task.list (required).
params.dir_in        = null

// Where Nextflow publishDir copies each task's results.
params.dir_out       = 'results'

// Shared files for every job. Default: ${dir_in}/template
params.template      = null

// File with one PALM argument line per task. Default: ${dir_in}/task.list
params.task_list     = null

// First task's -seed; task k uses (base_seed + k - 1). Must stay within
// PALM's accepted range: a positive integer up to 2^32.
params.base_seed     = 1

// Prefix used only in Nextflow job tags (for logs / monitoring), not in
// PALM output filenames (-o is set separately below).
params.report_prefix = 'palm'

// Apptainer/Singularity image name or path (see containers/Apptainer.def).
params.container     = 'palm.sif'

// Executable invoked inside the container. The image installs PALM under
// /opt/PALM; if you run without a container and `palm` is on PATH, pass
// --palm_bin palm
params.palm_bin      = '/opt/PALM/palm'

// If true, each task's output directory is packed to out_<k>.tar.gz and
// that archive is what gets published (instead of the directory itself).
params.pack_results  = false

//============================================================================
// Process: one independent PALM run
//============================================================================

process run_palm {
    // Shown in Nextflow logs / timeline; includes seed for auditability.
    tag { "${params.report_prefix}:${name}#${index}:seed${seed}" }

    // Resource defaults live under withLabel: 'palm_task' in nextflow.config
    label 'palm_task'

    // Hash inputs deeply so changed template file contents invalidate cache.
    cache 'deep'

    container params.container

    // Copy process outputs from the work directory into dir_out.
    publishDir params.dir_out, mode: 'copy'

    input:
    // All files from template/, staged into the work directory (basenames
    // visible to the shell / to relative paths in the task line).
    path template_files
    // index: 1-based task number
    // line:  PALM args from task.list (no -seed / -o)
    // seed:  value passed as -seed
    // name:  short label for the tag (from -i file, or task<k>)
    tuple val(index), val(line), val(seed), val(name)

    output:
    // Directory out_<k> or archive out_<k>.tar.gz, depending on pack_results.
    path("${params.pack_results ? "out_${index}.tar.gz" : "out_${index}"}"), emit: results

    script:
    def out_dir = "out_${index}"
    // Workflow owns -seed and -o so every job is unique and non-clobbering.
    // PALM writes its usual maps/CSVs under the -o prefix (out_<k>/palm_*).
    def palm_cmd = "${params.palm_bin} ${line} -seed ${seed} -o ${out_dir}/palm"
    if (params.pack_results) {
        """
        mkdir -p ${out_dir}
        ${palm_cmd}
        echo "${seed}" > ${out_dir}/seed.txt
        tar czf ${out_dir}.tar.gz ${out_dir}
        """
    } else {
        """
        mkdir -p ${out_dir}
        ${palm_cmd}
        echo "${seed}" > ${out_dir}/seed.txt
        """
    }
}

//============================================================================
// Workflow
//============================================================================

workflow {

    // -------------------------------------------------------------------------
    // Resolve paths
    // -------------------------------------------------------------------------
    if (!params.dir_in) {
        error 'Provide --dir_in (directory containing template/ and task.list)'
    }

    def template_dir = params.template ?: "${params.dir_in}/template"
    def task_list    = params.task_list ?: "${params.dir_in}/task.list"

    // -------------------------------------------------------------------------
    // Task channel
    // Each item: tuple(index, line, seed, name)
    //
    // We collect lines first, then assign indices, so seeds are a stable
    // function of line order in task.list (not of scheduling order).
    // -------------------------------------------------------------------------
    Channel
        .fromPath(task_list, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        // Skip blanks and shell-style comments.
        .filter { it && !it.startsWith('#') }
        .ifEmpty { error 'the task list is empty!' }
        .map { line ->
            // -seed / -o are reserved: injecting them twice would be ambiguous.
            if (line =~ /(^|\s)-seed(\s|=|$)/) {
                error "Task line must not include -seed (assigned by the workflow): ${line}"
            }
            if (line =~ /(^|\s)-o(\s|=|$)/) {
                error "Task line must not include -o (assigned by the workflow): ${line}"
            }
            line
        }
        .toList()
        .flatMap { lines ->
            lines.withIndex().collect { line, i ->
                def index = i + 1
                def seed  = params.base_seed + index - 1
                // PALM accepts a positive integer seed up to 2^32.
                if (seed < 1 || seed > 4294967295L) {
                    error "Seed ${seed} for task ${index} is outside 1..2^32-1; adjust --base_seed"
                }
                // Prefer a tag name from the first -i argument; else task<k>.
                def m = (line =~ /(?:^|\s)-i\s+(\S+)/)
                def name = m ? file(m[0][1]).simpleName : "task${index}"
                tuple(index, line, seed, name)
            }
        }
        .set { tasks_ch }

    // -------------------------------------------------------------------------
    // Template channel
    // collect() -> one list reused (broadcast) to every run_palm job.
    // -------------------------------------------------------------------------
    Channel
        .fromPath("${template_dir}/*", type: 'any', checkIfExists: true)
        .collect()
        .ifEmpty { error "the template folder is empty: ${template_dir}" }
        .set { template_ch }

    // -------------------------------------------------------------------------
    // Launch: for each task, stage template files and run PALM.
    //
    // Note: analyses that must share one PALM process (e.g. -corrmod / -corrcon
    // / NPC across modalities in a single call) belong on *one* task.list line,
    // not split across lines.
    // -------------------------------------------------------------------------
    run_palm(template_ch, tasks_ch)
}

// Modelines - must stay within the first or last 10 lines of the file
// kate: syntax Groovy;
// vim: syntax=groovy
// -*- mode: groovy;-*-
