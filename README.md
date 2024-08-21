In this repository, we illustrate with simulation that sex disparity in liver waitlist outcomes can persist, 
even if prioritization of waitlist candidates is based on a medical urgency score which by itself makes unbiased
predictions of pre-transplant mortality in males and females.

To this end, we generate realistic populations of liver transplantation candidates and donors. 
Pre-transplant survival is simulated to depend only on MELD. To allocate the donor livers to patients,
we prioritize candidates on MELD. Thus, we allocate solely on a medical urgency score which is not 
biased against females.

However, we assume that females are less likely to accept liver graft offers. We show that this results
in substantial sex disparities in waitlist outcomes, with females being both less likely to be transplanted than males,
and more likely to die on the waitlist. We then show that when a Cox proportional hazards model is fitted to the 
observed data with waitlist mortality as the outcome and female sex and MELD as adjustment variables,
female sex is insignificant. We explain why this is expected.

The conclusions of this are that (i) using a Cox proportional hazards model for waitlist mortality is a suitable approach
to develop unbiased scores for medical urgency, but (ii) scores developed in such a way (MELD 3.0 and GEMA-Na) cannot
correct for disparities which are the result of inequities in accessing transplantation.
