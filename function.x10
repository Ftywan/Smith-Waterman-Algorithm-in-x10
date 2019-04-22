def getBlockPosition(var i:Int, var j:Int):Rail[Int] {

	var position:Rail[Int];

	//divide the score matrix to blocks
	//get the num of blocks
	//size of the martix is (length+1n)(length2+1n)
	var numOfBlocksInHight = length1 / NUM_ROWS_IN_BLOCK + 1n;
	var numOfBlocksInWidth = length2 / NUM_COLS_IN_BLOCK + 1n;

	//get the final position
	var leftI:Int;
	var leftJ:Int;
	var rightI:Int;
	var rightJ:Int;

	//topleft position: never excel, need not consider the special case
	leftI = i * NUM_ROWS_IN_BLOCK + 1n;
	leftJ = j * NUM_COLS_IN_BLOCK + 1n;

	//special case: edge of the matrix
	if(i == numOfBlocksInHight-1n) {
		rightI = length1;
	}else {
		rightI = (i + 1n) * NUM_ROWS_IN_BLOCK + 1n;
	}

	if(j == numOfBlocksInWidth-1n) {
		rightJ = length2;
	}else {
		rightJ = (j + 1n) * NUM_COLS_IN_BLOCK + 1n;
	}

	position = [leftI, leftJ, rightI, rightJ];

	return position;
}