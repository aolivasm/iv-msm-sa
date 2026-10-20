# Data

This directory is reserved for local input data. Its contents are ignored by
Git to reduce the risk of publishing confidential or large files.

## Oregon Health Insurance Experiment application

The public-use Oregon Health Insurance Experiment data and replication files
are distributed by Harvard Dataverse under custom access terms:

<https://doi.org/10.7910/DVN/SJG1ED>

After accepting those terms, extract the downloaded archive as
`data/OHIE_Public_Use_Files/` and run:

```sh
Rscript scripts/prepare_ohie_application_data.R
```

The preparation script uses only the administrative emergency-department
outcome and its linked descriptive and Medicaid-enrollment files. Raw and
prepared person-level files must not be committed.

### Credit and citations

We acknowledge the Oregon Health Insurance Experiment investigators and the
Oregon Health Study Group for conducting the study and making public-use data
available. Analyses using these data should cite the dataset, as required by
its terms, and the relevant original studies:

- Finkelstein, Amy; Baicker, Katherine; Taubman, Sarah; Wright, Bill;
  Bernstein, Mira; Gruber, Jonathan; Allen, Heidi; Newhouse, Joseph P.;
  Schneider, Eric; and Zaslavsky, Alan (2018).
  *Replication Data for: Oregon Health Insurance Experiment*.
  Harvard Dataverse, V1. [doi:10.7910/DVN/SJG1ED](https://doi.org/10.7910/DVN/SJG1ED).
- Finkelstein, A., Taubman, S., Wright, B., Bernstein, M., Gruber, J.,
  Newhouse, J. P., Allen, H., Baicker, K., and the Oregon Health Study Group
  (2012). The Oregon Health Insurance Experiment: Evidence from the First Year.
  *The Quarterly Journal of Economics*, 127(3), 1057–1106.
  [doi:10.1093/qje/qjs020](https://doi.org/10.1093/qje/qjs020).
- Taubman, S. L., Allen, H. L., Wright, B. J., Baicker, K., and
  Finkelstein, A. N. (2014). Medicaid Increases Emergency-Department Use:
  Evidence from Oregon's Health Insurance Experiment.
  *Science*, 343(6168), 263–268.
  [doi:10.1126/science.1246183](https://doi.org/10.1126/science.1246183).

### Data-use conditions

Read and accept the provider's current terms on the linked dataset page before
downloading or using the data. Those terms require bibliographic citation and
prohibit commercial use, attempts to identify participants directly or
indirectly, and use of the materials to obtain information about or further
contact with participants known to the user. Applicable legal restrictions
also apply. This summary does not replace the provider's terms.

The repository's MIT license applies to its original software, not to the OHIE
data, and does not override the data-use conditions. Participant data are not
distributed here; obtain them from the provider. Keep raw and prepared records,
fitted application caches, and person-level scores out of public repositories.
The analyses in this repository are the responsibility of its authors and do
not imply endorsement by the OHIE investigators or data distributor.
