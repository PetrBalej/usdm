# usdm

## Uncertainty Analysis for Species Distribution Modelling

The usdm package provides a set of functions to support dealing with problematic situations in species distribution modelling (e.g., multicollinearity, positional uncertainty). 

To detect whether predictor variables are subjected to multicollinearity, you may use vif (variance inflation factor) metric, and some methods implemeted in this package including vifstep or vifcor (a stepwise procedure to identify collinear variables).

To detect whether positional uncertainty in species data may affect SDMs, different strategies are implemented through using either global or local spatial autocorrelation. You may check the following links for more information:


https://r-gis.net/?q=positional_uncertainty (using global spatial autocorrelation)

https://r-gis.net/?q=positional_uncertainty2 (using local spatial autocorrelation)



To develop species distribution models (SDMs), you may use the sdm package.

---

## Tracking Exclusion Chains in vifcor and vifstep (2.1-999)

When using `vifcor()` and `vifstep()` to identify and remove collinear predictors, the stepwise exclusion procedure is now fully tracked. Two new slots in the resulting `VIF` object provide complete audit trail and transitive chains:

### New VIF Object Slots

- **`exclusionLog`** (data.frame): A detailed record of each exclusion step with columns:
  - `step`: Order of exclusion (integer)
  - `excluded`: Name of the removed variable
  - `correlated_with`: Direct pairing partner (for `vifcor`) or most-correlated remaining variable (for `vifstep`)
  - `correlation`: Absolute correlation between the pair
  - `vif_excluded`: VIF of the excluded variable at that step
  - `vif_retained`: VIF of the retained variable (NA for `vifstep`)

- **`chains`** (list): Transitive exclusion chains for each final retained variable. A named list where each element represents a variable that "won" at least one pairing; the vector contains all excluded predictors attributed to it, ordered from most recently to most distally excluded.

### Example Usage

```r
file <- system.file("external/spain.tif", package="usdm")
r <- rast(file)
v1 <- vifcor(r, th=0.9)
```

The `show()` method displays:
```
---------- Exclusion log --------
 step excluded correlated_with correlation vif_excluded vif_retained
    1     Bio5           Bio10   0.9466061 2.558863e+12    1520.1380
    2    Bio10            Bio1   0.9288835 1.520137e+03     776.7311

---------- Exclusion chains (final variable <- removed predictors) --------
Bio1 <- Bio10 <- Bio5
```

### Accessing the Chains Programmatically

```r
# View the exclusion log
v1@exclusionLog

# Access the chains
v1@chains
# Output: list(Bio1 = c("Bio10", "Bio5"))
# Meaning: Bio1 is the final representative; 
#          it defeated Bio10, which had previously defeated Bio5
```

This enhancement enables full transparency and reproducibility in the variable selection process, essential for scientific reporting and model validation.
