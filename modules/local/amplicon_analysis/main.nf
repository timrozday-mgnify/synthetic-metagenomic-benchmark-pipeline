// Run the EBI-Metagenomics amplicon-analysis-pipeline (AAP) on the generated
// reads via a nested `nextflow run`. AAP is a full pipeline whose workflow is
// driven by params.input (no channel input), so we compose it as an ordinary
// process: write an AAP-format samplesheet and launch it once for a whole batch.
//
// Samples sharing a DB config (see the profile subworkflow's group key) are
// batched into ONE nested run — the AAP samplesheet is natively multi-row and
// AAP namespaces every output under aap_out/<sample_id>/, so we pay the nested
// Nextflow startup once per DB instead of once per sample.
//
// Runs on the host (executor 'local', no container) so it reuses the host
// nextflow + container engine — AAP manages its own containers/DBs internally.
// mapseq_databases and other AAP params come from the optional -c config.
// An optional custom PIMENTO primer library (params.aap_std_primer_library) is forwarded
// as --std_primer_library so PIMENTO matches against known primers instead of its bundled set.
process RUN_AAP {
    tag "${metas[0].database ?: 'community'} (${metas.size()})"
    label 'process_single'
    executor 'local'

    input:
    // Batched input: `metas` + `layout` describe every sample in the group; the DB
    // slots are group-invariant (same database ⇒ same files) so the subworkflow passes
    // one copy. Reads are NOT staged: layout carries their absolute host paths and this
    // process is executor 'local', so the nested AAP run reads them directly — this also
    // avoids sub_* basename collisions when many samples land in one task. -resume stays
    // correct because those paths embed the upstream workdir hash.
    // Optional path slots are stageAs'd to distinct fixed names so the shared NO_FILE
    // placeholder doesn't collide across slots. `use_built` selects the DB source.
    tuple val(metas), val(layout), val(use_built), path(aap_config, stageAs: 'aap_config_in'), path(mapseq_fasta, stageAs: 'mapseq_db.fasta'), path(mapseq_tax, stageAs: 'mapseq_db.tax'), path(mapseq_otu, stageAs: 'mapseq_db.otu'), path(mapseq_mscluster, stageAs: 'mapseq_db.mscluster'), val(rfam_cm), val(rfam_claninfo)

    output:
    // Glob (not the bare dir) so publishDir's saveAs sees each file as aap_out/<id>/...
    // and can route it to that sample's publish dir (see conf/modules.config).
    tuple val(metas), path("aap_out/**"), emit: results
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // DB config / engine configs / profile are group-invariant — read from metas[0].
    def meta       = metas[0]
    // A pipeline-built mapseq DB wins over the pass-through params.aap_config.
    def dbname     = meta.database ?: 'community'
    def db_cfg     = use_built ? '-c aap.config' : (params.aap_config ? "-c ${aap_config}" : '')
    // Engine + extra -c files come from meta (samplesheet aap_configs/aap_profile or the
    // params fallback, resolved to absolute paths in main.nf). DB config first, so later
    // files override earlier and the engine config wins. -profile only when requested.
    def extra_cfg  = (meta.aap_configs ?: []).collect { "-c ${file(it, checkIfExists: true)}" }.join(' ')
    def prof_arg   = meta.aap_profile ? "-profile ${meta.aap_profile}" : ''
    // Optional custom PIMENTO primer library. Global param, absolute host path (executor local,
    // nested run reads it directly). null => aap falls back to PIMENTO's bundled library.
    def primer_lib_arg = params.aap_std_primer_library ? "--std_primer_library ${file(params.aap_std_primer_library, type: 'dir', checkIfExists: true)}" : ''
    // One CSV row per sample: [id, single_end, fastq_1, fastq_2] (absolute read paths).
    // Emit one printf per row (leading indentation is harmless for commands, unlike a
    // heredoc body) so the samplesheet has no stray whitespace.
    assert layout.every { it[2] } : "RUN_AAP: empty fastq_1 in layout for ${metas*.id}"
    def sheet_cmds = layout.collect { id, se, fq1, fq2 ->
        "printf '%s,%s,%s,%s\\n' '${id}' '${fq1}' '${fq2}' '${se}' >> aap_samplesheet.csv"
    }.join('\n    ')
    """
    ${use_built ? "write_aap_config.py --name '${dbname}' --fasta ${mapseq_fasta} --tax ${mapseq_tax} --otu ${mapseq_otu} --mscluster ${mapseq_mscluster} --rfam-covariance-model '${rfam_cm}' --rfam-claninfo '${rfam_claninfo}' --output aap.config" : "true"}

    # AAP samplesheet: sample,fastq_1,fastq_2,single_end (one row per batched sample).
    printf 'sample,fastq_1,fastq_2,single_end\\n' > aap_samplesheet.csv
    ${sheet_cmds}
    # Fail loud if a row was dropped rather than silently profiling a subset.
    [ \$(( \$(grep -c . aap_samplesheet.csv) - 1 )) -eq ${layout.size()} ] || { echo "RUN_AAP: samplesheet row count != ${layout.size()}" >&2; exit 1; }

    nextflow run ebi-metagenomics/amplicon-analysis-pipeline \\
        -r ${params.aap_revision} \\
        ${prof_arg} \\
        --input aap_samplesheet.csv \\
        --outdir aap_out \\
        ${primer_lib_arg} \\
        ${db_cfg} ${extra_cfg}

    # Report contract is *.mseq.gz; AAP emits uncompressed .mseq. Idempotent (no-op if gz).
    find aap_out -name '*.mseq' -exec gzip -f {} +

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplicon-analysis-pipeline: ${params.aap_revision}
        nextflow: \$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    // Mirror AAP's per-sample namespacing (aap_out/<id>/...) for every batched sample.
    def stub_cmds = metas.collect { m ->
        "mkdir -p aap_out/${m.id}/taxonomy-summary && touch aap_out/${m.id}/taxonomy-summary/${m.id}.krona.txt"
    }.join('\n    ')
    """
    ${stub_cmds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplicon-analysis-pipeline: stub
        nextflow: stub
    END_VERSIONS
    """
}
