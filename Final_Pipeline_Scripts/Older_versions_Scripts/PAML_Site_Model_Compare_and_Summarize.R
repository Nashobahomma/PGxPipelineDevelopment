# Read the parsed PAML output table with the model 0, 1, 2, 7, 8, and 8a lnL
# values and perform a series of model comparisons to idetify the best fitting
# model for each CYP. Will write a new CSV with the existing PAML output values
# and additional columns: 
#	1) P-value for the test that yielded the best-fit model
# 	2) M2 vs M1 Chi-squared value (test statistic value)
#	3) M2 vs V1 P-value
#	4) M8 vs M7 Chi-squared value (test statistic value)
#	5) M8 vs M7 P-value
#	6) M8 vs M8a Chi-squared value (test statistic value)
#	7) M8 vs M8a P-value
#
# Takes three arguments:
#	1) Output of PAML parser Python script
#	2) Filename to save the full PAML table with model comparison results
#	3) Filename to save the "digest" PAML results (CYP name, gene-wide omega for each model, test statistcs and  P-values for each model comparison)
# 2023-06-23
# Amber Nashoba and Tom Kono

rm(list = ls())

# Take arguments
arguments <- commandArgs(TRUE)
paml_table_input <- arguments[1]
full_model_comp_out_filename <- arguments[2]
digest_model_comp_out_filename <- arguments[3]

# Read it as CSV
paml_table_data <- read.csv(paml_table_input, header=TRUE)

# Define a function for the M2 vs M1 comparison. Return the Chi-squared
# test statstic value and P-value for this comparison.
m2_vs_m1 <- function(cyp_gene)
{
	# Collect the lnL for model 1
	m1_lnl <- as.numeric(cyp_gene['Model1.lnL'])
	# Colelct lnL for model 2
	m2_lnl <- as.numeric(cyp_gene['Model2.lnL'])
	# Collect the num. parameters for model 1
	m1_np <- as.numeric(cyp_gene['Model1.np'])
	# Collect the num. parameters for model 2
	m2_np <- as.numeric(cyp_gene['Model2.np'])
	# Calculate the Chi-squared test statistic and the P-value
	#	For this comparison, it is -2 * (M1 - M2)
	ts <- -2 * (m1_lnl - m2_lnl)
	pval <- pchisq(ts, df=(m2_np-m1_np), lower.tail=FALSE)
	return(c(M2.vs.M1.TestStatistic=ts, M2.vs.M1.Pvalue=pval))
}


# Define a function for M8 vs M7
m8_vs_m7 <- function(cyp_gene)
{
	m8_lnl <- as.numeric(cyp_gene['Model8.lnL'])
	m7_lnl <- as.numeric(cyp_gene['Model7.lnL'])
	m8_np <- as.numeric(cyp_gene['Model8.np'])
	m7_np <- as.numeric(cyp_gene['Model7.np'])
	# Calculate the Chi-squared test statistic and the P-value
	#	For this comparison, it is -2 * (M7 - M8)
	ts <- -2 * (m7_lnl - m8_lnl)
	pval <- pchisq(ts, df=(m8_np-m7_np), lower.tail=FALSE)
	return(c(M8.vs.M7.TestStatistic=ts, M8.vs.M7.Pvalue=pval))
}

# Define a function for M8 vs M8a
m8_vs_m8a <- function(cyp_gene)
{
	m8_lnl <- as.numeric(cyp_gene['Model8.lnL'])
	m8a_lnl <- as.numeric(cyp_gene['Model8a.lnL'])
	m8_np <- as.numeric(cyp_gene['Model8.np'])
	m8a_np <- as.numeric(cyp_gene['Model8a.np'])
	# Calculate the Chi-squared test statistic and the P-value
	#	For this comparison, it is -2 * (M8a - M8)
	ts <- -2 * (m8a_lnl - m8_lnl)
	pval <- pchisq(ts, df=(m8_np-m8a_np), lower.tail=FALSE)
	return(c(M8.vs.M8a.TestStatistic=ts, M8.vs.M8a.Pvalue=pval))
}

# Apply the three model comparison functions over the whole table
m2_m1_comp <- apply(paml_table_data, 1, m2_vs_m1)
m8_m7_comp <- apply(paml_table_data, 1, m8_vs_m7)
m8_m8a_comp <- apply(paml_table_data, 1, m8_vs_m8a)

# "cbind()" the results onto the original table. Have to transpose the
# outputs of the model test functions 
paml_w_model_comparisons <- cbind(paml_table_data, t(m2_m1_comp))
paml_w_model_comparisons <- cbind(paml_w_model_comparisons, t(m8_m7_comp))
paml_w_model_comparisons <- cbind(paml_w_model_comparisons, t(m8_m8a_comp))

# Save the full output table
write.csv(paml_w_model_comparisons, file=full_model_comp_out_filename, row.names=FALSE, quote=TRUE)

# Save the "digest"
digested_paml_table <- paml_w_model_comparisons[, 
	c("CYP.Name", 
	"Model1.dNdS", "Model2.dNdS", "M2.vs.M1.TestStatistic",  "M2.vs.M1.Pvalue",
	"Model7.dNdS", "Model8.dNdS", "M8.vs.M7.TestStatistic",  "M8.vs.M7.Pvalue",
	"Model8a.dNdS", "M8.vs.M8a.TestStatistic",  "M8.vs.M8a.Pvalue")]
write.csv(digested_paml_table, file=digest_model_comp_out_filename, row.names=FALSE, quote=TRUE)