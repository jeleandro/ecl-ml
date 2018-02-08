IMPORT ML;
/* *********************************************************************************************************************

Mahalanobis.ecl - new implementation
-fixes:
-- Lower Gamma function was replaced
-- Standartization was using populational standard deviation and was changed to sample std.


**********************************************************************************************************************/
ChiSquareCDF(INTEGER df, REAL cv) := FUNCTION
	RETURN ML.Utils.lowerGamma(df / 2, cv / 2) / ML.Utils.gamma(df / 2);
END;


Export Mahalanobis(DATASET(ML.Mat.Types.Element) A, REAL sensitivity = 0.05) := MODULE
	ZComp := ML.Mat.Pca(A).ZComp;

	// Calculate sample standard deviation for each column
	stdev := TABLE(ZComp, {
		y;
		stdev:= SQRT(	sum(GROUP, value * value)/(max(GROUP, x)-1 ));
	}, y);

	// Creates ZComp with normalized columns
	ZCompNorm := PROJECT(JOIN(ZComp, stdev, LEFT.y = RIGHT.y), TRANSFORM(
		ML.Mat.Types.Element,
		SELF.x := LEFT.x;
		SELF.y := LEFT.y;
		SELF.value := LEFT.value / LEFT.stdev;
	));

	// Calculates squared distance
	EXPORT dsq := PROJECT(
		TABLE(ZCompNorm, {x;dsq := SUM(GROUP, value * value);}, x),
		TRANSFORM(
			ML.Mat.Types.Element,
			SELF.x := LEFT.x;
			SELF.y := 1;
			SELF.value := LEFT.dsq;
		)
	);

	// Calculates how many degrees of freedom we are dealing with
	df := ML.Mat.Has(A).Stats.YMax;
	

	// turn distances into probabilities
	EXPORT prob := PROJECT(dsq, TRANSFORM(
		ML.Mat.Types.Element,
		SELF.x := LEFT.x;
		SELF.y := 1;
		SELF.value := ChiSquareCDF(df, LEFT.value);
	));	
	

	// Verifies if p-value of each DSQ exceeds defined sensitivity
	EXPORT is_outlier := PROJECT(prob, TRANSFORM(
		ML.Mat.Types.Element,
		SELF.x := LEFT.x;
		SELF.y := 1;
		SELF.value := IF(LEFT.value > (1 - sensitivity), 1, 0);
	));
END;