After phrase *Insert method of multiple correction*, program expects one of the following methods:

<li> bonferroni : one-step correction </li>

<li> sidak : one-step correction </li>

<li> holm-sidak : step down method using Sidak adjustments </li>

<li> holm : step-down method using Bonferroni adjustments </li>

<li> simes-hochberg : step-up method (independent) </li>

<li> hommel : closed method based on Simes tests (non-negative) </li>

<li> fdr_bh : Benjamini/Hochberg (non-negative) </li>

<li> fdr_by : Benjamini/Yekutieli (negative) </li>

<li> fdr_tsbh : two stage fdr correction (non-negative) </li>

<li> fdr_tsbky : two stage fdr correction (non-negative) </li>

<br>
If input is empty, returns p-values without correction.
