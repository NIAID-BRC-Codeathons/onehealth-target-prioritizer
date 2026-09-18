# OneHealth Conserved Target Prioritizer

**NIAID-BRCs AI Codeathon 2.0** · September 16–18, 2026 · Argonne National Laboratory

Identifying conserved vaccine or therapeutic targets for an emerging or re-emerging virus using genomic surveillance, immune evidence, protein structure, host interaction, and literature data.

Project page: https://niaid-brc-codeathons.github.io/projects/onehealth-target-prioritizer/

---

> **This is a draft pitch, not a plan.**
>
> What follows is a one-slide proposal from the organizing team. It exists
> to seed a team, not to constrain one. Scope, methods, target organism,
> and success criteria are all still open — expect them to change
> substantially. Turning this into a real plan is the team's first job, and
> it lands in the project charter due August 28, 2026.

---

## Goal (proposed)

Identify conserved vaccine or therapeutic targets for an emerging or re-emerging virus using genomic surveillance, immune evidence, protein structure, host interaction, and literature data.

## Three-Day MVP (proposed)

Select one virus family and one narrowly defined product objective. Retrieve human and animal sequences, measure conservation and lineage coverage, annotate epitopes and functional domains, identify host interactions and escape mutations, and rank 10–20 candidate targets.

The final output should be a transparent target dossier rather than a claim that a vaccine or antiviral has been designed.

## Evaluation (proposed)

Recovery of established targets, conservation across host species and lineages, evidence density, structural accessibility, and stability of rankings under alternative weighting schemes.

## Leads

- Maged Hemida
- Rong Hai

Team assignments are still being finalized. Participants can review their project, and request a reassignment, in the participant spreadsheet circulated by the organizing team.

## Working here

This repository is the team's working space for the codeathon — code, notebooks, data pointers, and notes. Replace this README with the real thing once the charter is written. Team members get access through the [NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization; accept the invitation if you have not already.

## In this repo

- [`integration_retrieve_qc_align_conserve_match/`](integration_retrieve_qc_align_conserve_match/) — end-to-end viral protein analysis workflow that connects NCBI sequence retrieval, JSON-to-FASTA conversion, MAFFT multiple-sequence alignment, pairwise identity analysis, conservation analysis, conserved peptide-region extraction, and peptide matching against `cov_pos_ab.tsv`. The integrated runner accepts a CSV containing NCBI `Taxonomy ID` values, creates run-specific outputs under `pipeline_runs/`, and coordinates the component scripts from retrieval through peptide matching. See [`integration_retrieve_qc_align_conserve_match/README.md`](integration_retrieve_qc_align_conserve_match/README.md) for dependencies, setup, command-line options, workflow stages, output files, and troubleshooting.

- [`conservation/`](conservation/) — standalone CLI tool that finds conserved
  peptide-length windows in a protein multiple sequence alignment (input: an
  existing MSA; this repo does not align sequences). Use this component when
  you already have an MSA and only need conserved-region analysis. See
  [`conservation/README.md`](conservation/README.md) for usage and
  [`conservation/PIPELINE_INTEGRATION.md`](conservation/PIPELINE_INTEGRATION.md)
  for its input/output contract.
