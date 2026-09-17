# JSON Protein Records to FASTA Converter

A small Python command-line utility for converting protein sequence
records from a text file containing **multiple concatenated JSON
objects** into a standard FASTA file.

The input format is useful for data exports where each JSON object
contains metadata such as a taxonomic ID and a `records` array. Each
protein record is expected to contain fields such as `accession`,
`title`, `organism`, `length`, and `sequence`.

## Why this script is needed

The input file is not necessarily a single valid JSON document. It may
contain several complete JSON objects one after another, for example:

``` json
{
  "taxid": 2697049,
  "records": [
    {
      "accession": "YP_009742617.1",
      "title": "nsp10 [Severe acute respiratory syndrome coronavirus 2]",
      "length": 139,
      "organism": "Severe acute respiratory syndrome coronavirus 2",
      "sequence": "AGNATEVPANSTVLSFCAFAVDAAKAY..."
    }
  ]
}
{
  "taxid": 694069,
  "records": []
}
```

Because this is not one JSON array or object, using `json.load()`
directly produces an error similar to:

``` text
json.decoder.JSONDecodeError: Extra data
```

The converter avoids this problem by using
`json.JSONDecoder().raw_decode()` repeatedly to parse each JSON object
separately.

## Requirements

-   Python 3.7 or later
-   No third-party Python packages are required

The script uses only modules from the Python standard library:

-   `argparse`
-   `json`
-   `re`
-   `pathlib`

## Usage

Assuming the converter is saved as `json_to_fasta.py` and the input file
is `out_complete.txt`:

``` bash
python json_to_fasta.py out_complete.txt -o converted.fasta
```

The general syntax is:

``` bash
python json_to_fasta.py INPUT_FILE [-o OUTPUT_FILE] [-w WIDTH]
```

### Arguments

  -----------------------------------------------------------------------
  Argument                            Description
  ----------------------------------- -----------------------------------
  `input`                             Input text file containing one or
                                      more concatenated JSON objects.

  `-o`, `--output`                    Output FASTA filename. Default:
                                      `converted.fasta`.

  `-w`, `--width`                     Number of amino-acid residues per
                                      FASTA sequence line. Default: `60`.
  -----------------------------------------------------------------------

For example, to wrap sequences at 80 residues per line:

``` bash
python json_to_fasta.py out_complete.txt -o proteins.fasta --width 80
```

## Expected input record

Each sequence record should contain at least an accession and sequence:

``` json
{
  "uid": "1826688927",
  "accession": "YP_009742617.1",
  "title": "nsp10 [Severe acute respiratory syndrome coronavirus 2]",
  "length": 139,
  "organism": "Severe acute respiratory syndrome coronavirus 2",
  "sequence": "AGNATEVPANSTVLSFCAFAVDAAKAY..."
}
```

The fields used by the converter are:

-   `accession` --- required; used as the primary FASTA identifier.
-   `sequence` --- required; written as the FASTA sequence.
-   `title` --- optional; included in the FASTA header.
-   `organism` --- optional; included in the FASTA header.
-   `length` --- optional; compared with the actual sequence length for
    validation.

Other fields, such as `uid`, are ignored.

## FASTA output

A record is written in the following form:

``` text
>YP_009742617.1 nsp10 [Severe acute respiratory syndrome coronavirus 2] | Severe acute respiratory syndrome coronavirus 2
AGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPE
ANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVC
GMWKGYGCSCDQLREPMLQ
```

By default, sequences are wrapped at 60 residues per line.

## Validation and error handling

The converter performs several checks while processing the input.

### Missing accession or sequence

Records without an `accession` or `sequence` are skipped and reported to
the terminal.

### Sequence length

If a record contains a `length` field, the script compares it with the
number of residues in the cleaned sequence. A mismatch produces a
warning such as:

``` text
WARNING: ABC123.1: reported length=200, actual length=198
```

The record is still written to the FASTA file.

### Invalid JSON

If the script encounters content that cannot be decoded as a JSON
object, it reports the approximate character position, line, and column
where parsing failed.

### Empty `records` arrays

JSON objects containing an empty `records` array are valid and simply
contribute no FASTA entries.

## Sequence cleaning

Before writing a sequence, the script:

1.  removes whitespace;
2.  converts all residues to uppercase; and
3.  wraps the resulting sequence according to the requested line width.

The script does not otherwise alter, translate, or filter amino-acid
characters.

## Duplicate records

The converter intentionally does **not** remove duplicate accessions or
identical sequences. Every valid record encountered in the input is
written to the FASTA output.

This behavior preserves the source data and avoids silently discarding
records. Deduplication, if required for downstream analysis, should be
performed as a separate step using explicitly chosen criteria.

## Example terminal output

A successful run may look like:

``` text
Conversion complete
-------------------
JSON objects read : 10
FASTA records     : 429
Records skipped   : 0
Length warnings   : 0
Output            : converted.fasta
```

The exact counts depend on the input file.

## Troubleshooting

### `Extra data` error

If you see:

``` text
Invalid JSON: Extra data
```

you are likely running an earlier version of the converter that uses
`json.load()`. Use the version that parses concatenated JSON objects
with `json.JSONDecoder().raw_decode()`.

### `Input file not found`

Check that the input filename and path are correct:

``` bash
python json_to_fasta.py /full/path/to/out_complete.txt -o converted.fasta
```

### No sequence records found

Confirm that the JSON objects contain a `records` field and that it is
an array of protein records.

## Reproducible example

Keep the script and source data in the same directory:

``` text
project/
├── json_to_fasta.py
├── out_complete.txt
└── README.md
```

Then run:

``` bash
cd project
python json_to_fasta.py out_complete.txt -o converted.fasta
```

After conversion:

``` text
project/
├── json_to_fasta.py
├── out_complete.txt
├── converted.fasta
└── README.md
```

## Notes

-   The converter is intended for protein sequence records but the
    FASTA-writing logic is generic.
-   The input may contain any number of concatenated JSON objects.
-   JSON objects without sequence records are safely ignored.
-   The original input file is not modified.
-   README generated and script co-authored by by GPT-5.6 SoI
