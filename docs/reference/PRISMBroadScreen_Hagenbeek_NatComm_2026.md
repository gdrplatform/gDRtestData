# PRISM Broad Screen Cancer Cell Line Viability Dataset

The single-agent screening data from the PRISM platform - with one drug
and 774 cell lines.

## Format

A `MultiAssayExperiment` object containing a single-agent experiment
with 5 assays.

## Source

Hagenbeek et al, Nat Comm 2026, Accepted

## Loading

`gDRutils::get_synthetic_data("prism")`

## Description

### Assays

Experimental data matrices tracking cell line viability across
processing stages:

- RawTreated:

  Raw intensity values for treated cell line pools.

- Controls:

  Raw intensity values for DMSO and killing controls.

- Normalized:

  Viability scores relative to controls.

- Averaged:

  Replicate-averaged viability scores.

- Metrics:

  Dose-response curve parameters (e.g., IC50, AUC).

### Column Metadata (colData)

Feature metadata describing the cell lines (models) across assays:

- clid:

  Unique Broad Institute cell line identifier.

- CellLineName:

  Publicly recognized cancer cell line name.

- Tissue:

  Primary tissue of origin.

- parental_identifier:

  Identifier for the parental cell line.

- subtype:

  Histological or molecular subtype.

- ReferenceDivisionTime:

  Estimated doubling time in hours.

### Row Metadata (rowData)

Feature metadata describing the treatment compounds across assays:

- Gnumber:

  Genentech compound identifier.

- DrugName:

  Public or commercial name of the tested compound.

- drug_moa:

  Mechanism of action / biological target of the drug.

- Duration:

  Exposure time of the cell lines to the compound (e.g., 72h).
