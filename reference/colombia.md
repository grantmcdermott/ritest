# Vendor data from an RCT conducted in Bogota, Colombia.

This dataset is a preliminary version of the data described in Iacovone
and McKenzie (forthcoming). It conforms to an earlier blog post written
by one of the authors (McKenzie, 2017) and contains data collected
during a randomized controlled trial on supply chains among fresh
produce vendors in Bogota, Colombia.

## Usage

``` r
colombia
```

## Format

A data frame with 3336 rows and 8 variables:

- b_block:

  business block (geographic identifier)

- b_pair:

  randomly assigned treatment pairs (plus one triplet)

- b_treat:

  treatment status (1 = yes, 2 = no)

- dayscorab:

  number of days that week requiring visits to Corabastos central market

- b_dayscorab:

  baseline of the previous variable

- miss_b_dayscorab:

  dummy for missing baseline information (1 = yes, 2 = no)

- round2, round3:

  survey round dummies

## Details

The RCT studies the impact of purchase aggregation by many
microenterprises (here: fruit and vegetable vendors), which enables a
reduction in costly individual visits to a large central market. My
thanks to David McKenzie for sharing the data.

## References

McKenzie, D. (2017) "Finally, a way to do easy randomization inference
in Stata!", Development Impact (World Bank blog).
<https://blogs.worldbank.org/impactevaluations/finally-way-do-easy-randomization-inference-stata>.

Iacovone, L. and McKenzie, D. (forthcoming) "Shortening Supply Chains:
Experimental Evidence from Fruit and Vegetable Vendors in Bogota",
Economic Development and Cultural Change.
<https://doi.org/10.1086/714050>.
