
function rand_matrix(rows, cols) {
    var C = new Array(rows);  
    for (var r = 0; r < rows; ++r) {
	C[r] = new Array(cols); 
	for (var c = 0; c < cols; ++c) {
	    C[r][c] = Math.random();
	}
    }
    return C;
}

function create_matrix(rows, cols) {
    var C = new Array(rows);  
    for (var r = 0; r < rows; ++r) {
	C[r] = new Array(cols); 
	for (var c = 0; c < cols; ++c) {
	    C[r][c] = 0;
	}
    }
    return C;
}

function mult(C, A, B) {
    var aNumRows = A.length
    var aNumCols = A[0].length
    var bNumRows = B.length
    var bNumCols = B[0].length

    for (var r = 0; r < aNumRows; ++r) {
	for (var c = 0; c < bNumCols; ++c) {
	    C[r][c] = 0;             
	    for (var i = 0; i < aNumCols; ++i) {
		C[r][c] += A[r][i] * B[i][c];
	    }
	}
    }

    return C;
}

exports.rand_matrix= rand_matrix;
exports.create_matrix= create_matrix;
exports.mult= mult;
