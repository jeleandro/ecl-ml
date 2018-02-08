IMPORT ML;

/* ********************************************************************************************************************* */

/**
 * This function ranks all elments across the cells of the matrix.
 * WARNING: filter the data to rank only one variable
 * Param: A => the matrix of elements
 * Param: groups => desired number of groups if none was especified the it is assumed the number of elements.
 * Param: ties   => in case of 2 or more element have the same value ties specify how to deal with that.
 *                  The default is do nothing and assumed different ranks. It is also possible l=mininum rank, m=mean of the ranks, h=maximum
 *
 */
export RankElements(DATASET(ML.Mat.Types.Element) A, UNSIGNED groups = 0,STRING1 ties='l') := FUNCTION
	sortedA := SORT(A, value);
	N := ML.Mat.Has(sortedA).Stats.NElements;
	fgroups := IF(groups = 0, N, groups);
	
	dense := PROJECT(
		sortedA,
		TRANSFORM(
		  {ML.Mat.Types.Element; real rankOrder},
			SELF.rankOrder := TRUNCATE(COUNTER * 1000 / (N + 1)),
			SELF := LEFT
		)
	);
	
	tiesGroups := TABLE(dense,{
			value,
			rankLow  := min(group,rankOrder),
			rankMean := ave(group,rankOrder),
			rankHigh := max(group,rankOrder)
		},
		value
	);  
	ranks := JOIN (dense, tiesGroups, LEFT.value = RIGHT.value, TRANSFORM(
		  ML.Mat.Types.Element,
		  self.x := LEFT.x;
		  self.y := LEFT.y;
		  self.value := map(
					ties = 'l' => RIGHT.rankLow,
					ties = 'm' => RIGHT.rankMean,
					ties = 'h' => RIGHT.rankHigh,
					ties = 'o' => LEFT.rankOrder,
					LEFT.rankOrder
			)
		)
	);
	return ranks;
END;