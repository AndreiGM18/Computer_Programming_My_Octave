// Mitran Andrei-Gabriel, 313CA
#ifndef MY_OCTAVE_FUNC
#define MY_OCTAVE_FUNC

// This defines the value which helps prevent stack overflow.
#define MOD 10007

// The following functions print an error message related to memory
// allocation.
void error_alloc1(int *p)
{
	if (!p) {
		fprintf(stderr, "Could not allocate memory\n");
		exit(-1);
	} else {
		return;
	}
}

void error_alloc2(int **p)
{
	if (!p) {
		fprintf(stderr, "Could not allocate memory\n");
		exit(-1);
	} else {
		return;
	}
}

void error_alloc3(int ***p)
{
	if (!p) {
		fprintf(stderr, "Could not allocate memory\n");
		exit(-1);
	} else {
		return;
	}
}

// This function checks if a matrix with the given index exists.
// It also prints an error message if necessary.
int check_index(int i, int count_matrix)
{
	if (i >= count_matrix || i < 0) {
		printf("No matrix with the given index\n");
		return 0;
	} else {
		return 1;
	}
}

// This function allocates more memory to the vector of matrixes
// and to the arrays which hold the number of rows and columns.
void alloc_vector(int ****matrix, int **m, int **n, int count_matrix)
{
	// Allocating more memory to the vector of matrixes
	int ***temp_matrix = realloc(*matrix, (count_matrix + 1) * sizeof(int **));
	error_alloc3(temp_matrix);

	// Allocating more memory to the array
	// which holds the number of rows of each matrix
	int *temp_m = realloc(*m, (count_matrix + 1) * sizeof(int));
	error_alloc1(temp_m);

	// Allocating more memory to the array
	// which holds the number of columns of each matrix
	int *temp_n = realloc(*n, (count_matrix + 1) * sizeof(int));
	error_alloc1(temp_n);

	// Tying the original ones to the temporary ones
	*matrix = temp_matrix;
	*m = temp_m;
	*n = temp_n;
}

// This function gets a matrix from stdin.
int **read(int row, int column)
{
	// Creating the pointer to the new matrix
	int **mat;

	// Allocating memory for the matrix
	mat = malloc(row * sizeof(int *));
	error_alloc2(mat);

	// Allocating memory for each line
	for (int i = 0; i < row; ++i) {
		mat[i] = malloc(column * sizeof(int));
		error_alloc1(mat[i]);
	}

	// Reading each element
	for (int i = 0; i < row; ++i)
		for (int j = 0; j < column; ++j)
			scanf("%d", &mat[i][j]);

	// Returning the pointer
	return mat;
}

// This function adds a matrix to the vector.
void add_to_vector(int ****matrix, int **m, int **n, int *count_matrix)
{
	// Declaring and reading the number of rows and columns
	int row, column;
	scanf("%d %d", &row, &column);

	// Calling this function to allocate more memory to the vector
	alloc_vector(matrix, m, n, *count_matrix);

	// Reading the new matrix and adding it to the vector
	(*matrix)[*count_matrix] = read(row, column);

	// Adding the number of rows and columns to the arrays that hold
	// the two dimensions
	(*m)[*count_matrix] = row;
	(*n)[*count_matrix] = column;

	// Adding one to the number of matrixes in the vector
	++(*count_matrix);
}

// This function prints the dimensions of a specified matrix.
void show_dim(int *m, int *n, int count_matrix)
{
	// Declaring its index and reading it
	int index;
	scanf("%d", &index);

	// Checking it and printing the dimensions if possible
	if (check_index(index, count_matrix))
		printf("%d %d\n", m[index], n[index]);
}

// This function prints a specified matrix.
void show(int ***matrix, int *m, int *n, int count_matrix)
{
	// Declaring its index and reading it
	int index;
	scanf("%d", &index);

	// Checking the index and printing the matrix if possible
	if (check_index(index, count_matrix)) {
		for (int i = 0; i < m[index]; ++i) {
			for (int j = 0; j < n[index]; ++j)
				printf("%d ", matrix[index][i][j]);
			printf("\n");
		}
	}
}

// This function resizes a specified matrix.
void resize(int ****matrix, int **m, int **n, int count_matrix)
{
	// Declaring its index and reading it
	int index;
	scanf("%d", &index);

	// Checking the index
	if (check_index(index, count_matrix)) {
		// Declaring the new number of rows and columns and two temporary
		// arrays, which will hold the indexes of each element that will
		// appear in the resized matrix
		int new_row, new_col;
		int *r, *c, **mat;

		// Reading the new number of rows
		scanf("%d", &new_row);

		// Allocating memory and reading the row index of each new element
		r = malloc(new_row * sizeof(int));
		error_alloc1(r);
		for (int i = 0; i < new_row; ++i)
			scanf("%d", &r[i]);

		// Reading the new number of columns
		scanf("%d", &new_col);

		// Allocating memory and reading the column index of each new element
		c = malloc(new_col * sizeof(int));
		error_alloc1(c);
		for (int i = 0; i < new_col; ++i)
			scanf("%d", &c[i]);

		// Allocating memory for the resized matrix
		mat = malloc(new_row * sizeof(int *));
		error_alloc2(mat);
		for (int i = 0; i < new_row; ++i) {
			mat[i] = malloc(new_col * sizeof(int));
			error_alloc1(mat[i]);
		}

		// Adding the specified element to the resized matrix
		for (int i = 0; i < new_row; ++i)
			for (int j = 0; j < new_col; ++j)
				mat[i][j] = (*matrix)[index][r[i]][c[j]];

		// Freeing the memory allocated for the original matrix
		for (int i = 0; i < (*m)[index]; ++i)
			free((*matrix)[index][i]);
		free((*matrix)[index]);

		// Adding the resized matrix in the original's place
		(*matrix)[index] = mat;

		// Changing the number of rows and columns in the arrays
		// that hold them
		(*m)[index] = new_row;
		(*n)[index] = new_col;

		// Freeing the temporary arrays
		free(r);
		free(c);
	}
}

// This functions adds two matrixes.
void add(int **A, int **B, int ***C, int dim)
{
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			(*C)[i][j] = (A[i][j] + B[i][j] + MOD) % MOD;
}

// This function subtracts the second matrix from the first.
void sub(int **A, int **B, int ***C, int dim)
{
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			(*C)[i][j] = (A[i][j] - B[i][j] + 2 * MOD) % MOD;
}

// This function allocates the required memory for one of the smaller
// matrixes that are the result of partitioning the bigger one.
void alloc_new_dim(int ***M, int new_dim)
{
	int **temp_M;
	temp_M = malloc(new_dim * sizeof(int *));
	error_alloc2(temp_M);
	for (int i = 0; i < new_dim; ++i) {
		temp_M[i] = malloc(new_dim * sizeof(int));
		error_alloc1(temp_M[i]);
	}
	*M = temp_M;
}

// This function allocates memory for some variables that are required for
// the Strassen algorithm when a recursive call is made.
void alloc_all(int ****a, int ****b, int ****c, int ****M, int new_dim)
{
	int ***temp_a = malloc(4 * sizeof(int **));
	error_alloc3(temp_a);
	int ***temp_b = malloc(4 * sizeof(int **));
	error_alloc3(temp_b);
	int ***temp_c = malloc(4 * sizeof(int **));
	error_alloc3(temp_c);
	int ***temp_M = malloc(7 * sizeof(int **));
	error_alloc3(temp_M);
	for (int i = 0; i < 4; ++i) {
		alloc_new_dim(&temp_a[i], new_dim);
		alloc_new_dim(&temp_b[i], new_dim);
		alloc_new_dim(&temp_c[i], new_dim);
		}
	for (int i = 0; i < 7; ++i)
		alloc_new_dim(&temp_M[i], new_dim);
	*a = temp_a;
	*b = temp_b;
	*c = temp_c;
	*M = temp_M;
}

// This function frees the allocated memory for one of the smaller
// matrixes that are the result of partitioning a bigger one.
void free_new_dim(int **M, int new_dim)
{
	for (int i = 0; i < new_dim; ++i)
		free(M[i]);
	free(M);
}

// This function frees the allocated memory for the variables that were
// required for the Strassen algorithm when the recursive call was made.
void free_all(int ***a, int ***b, int ***c, int ***M, int new_dim)
{
	for (int i = 0; i < 4; ++i) {
		free_new_dim(a[i], new_dim);
		free_new_dim(b[i], new_dim);
		free_new_dim(c[i], new_dim);
		}
	free(a);
	free(b);
	free(c);
	for (int i = 0; i < 7; ++i)
		free_new_dim(M[i], new_dim);
	free(M);
}

// This function partitions a bigger matrix into four smaller ones.
void part(int ****m, int **M, int new_dim)
{
	for (int i = 0; i < new_dim; ++i)
		for (int j = 0; j < new_dim; ++j) {
			(*m)[0][i][j] = M[i][j];
			(*m)[1][i][j] = M[i][j + new_dim];
			(*m)[2][i][j] = M[i + new_dim][j];
			(*m)[3][i][j] = M[i + new_dim][j + new_dim];
		}
}

// This function constructs a matrix from four smaller ones.
void construct(int ***m, int ***M, int new_dim)
{
	for (int i = 0; i < new_dim; i++) {
		for (int j = 0; j < new_dim; j++) {
			(*M)[i][j] = m[0][i][j];
			(*M)[i][j + new_dim] = m[1][i][j];
			(*M)[i + new_dim][j] = m[2][i][j];
			(*M)[i + new_dim][j + new_dim] = m[3][i][j];
			}
		}
}

// This function represents the Strassen algorithm. It makes
// a recursive call to itself every time a multiplication
// is required. The two matrixes that are to be multiplied
// are called A and B. The result is C. The partitions
// for any of them have 11, 12, 21 or 21 after their name,
// meaning upper left, upper right, lower left or lower right
// respectively. Since the matrixes are guaranteed to be 2^n x 2^n,
// let dim be 2^n.
void strassen_algorithm(int **A, int **B, int ***C, int dim)
{
	if (dim == 1) {
		// Multiplying, if after partitioning
		// the resulting matrixes are 1 x 1
		(*C)[0][0] = (A[0][0] * B[0][0] + MOD * MOD) % MOD;
	} else {
		// Dividing the current dimension by 2
		// Declaring the required variables
		int new_dim = dim / 2;
		int ***a, ***b, ***c, ***M;
		int **result_1, **result_2;

		// Allocating memory
		alloc_all(&a, &b, &c, &M, new_dim);
		alloc_new_dim(&result_1, new_dim);
		alloc_new_dim(&result_2, new_dim);

		// Partitioning A and B
		part(&a, A, new_dim);
		part(&b, B, new_dim);

		// Calclulating M1 = (A11 + A21) * (B11 + B21)
		add(a[0], a[3], &result_1, new_dim);
		add(b[0], b[3], &result_2, new_dim);
		strassen_algorithm(result_1, result_2, &M[0], new_dim);

		// Calculating M2 = (A21 + A22) * B11
		add(a[2], a[3], &result_1, new_dim);
		strassen_algorithm(result_1, b[0], &M[1], new_dim);

		// Calculating M3 = A11 * (B12 - B22)
		sub(b[1], b[3], &result_2, new_dim);
		strassen_algorithm(a[0], result_2, &M[2], new_dim);

		// Calculating M4 = A22 * (B21 - B11)
		sub(b[2], b[0], &result_2, new_dim);
		strassen_algorithm(a[3], result_2, &M[3], new_dim);

		// Calculating M5 = (A11 + A12) * B22
		add(a[0], a[1], &result_1, new_dim);
		strassen_algorithm(result_1, b[3], &M[4], new_dim);

		// Calculating M6 = (A21 - A11) * (B11 + B12)
		sub(a[2], a[0], &result_1, new_dim);
		add(b[0], b[1], &result_2, new_dim);
		strassen_algorithm(result_1, result_2, &M[5], new_dim);

		// Calculating M7 = (A12 - A22) * (B21 + B22)
		sub(a[1], a[3], &result_1, new_dim);
		add(b[2], b[3], &result_2, new_dim);
		strassen_algorithm(result_1, result_2, &M[6], new_dim);

		// Calculating C12 = M3 + M5
		add(M[2], M[4], &c[1], new_dim);

		// Calculating C21 = M2 + M4
		add(M[1], M[3], &c[2], new_dim);

		// Calculating C11 = M1 + M4 - M5 + M7
		add(M[0], M[3], &result_1, new_dim);
		add(result_1, M[6], &result_2, new_dim);
		sub(result_2, M[4], &c[0], new_dim);

		// Calculating C22 = M1 - M2 + M3 + M6
		add(M[0], M[2], &result_1, new_dim);
		add(result_1, M[5], &result_2, new_dim);
		sub(result_2, M[1], &c[3], new_dim);

		// Grouping the results obtained in a single matrix
		construct(c, C, new_dim);

		// Freeing the allocated memory
		free_all(a, b, c, M, new_dim);
		free_new_dim(result_1, new_dim);
		free_new_dim(result_2, new_dim);
	}
}

// This function utilizes the Strassen multiplication algorithm and returns
// a pointer to the result.
int **strassen_multiply(int ****matrix, int **m, int i1, int i2)
{
	int dim = (*m)[i1];
	int **C;
	alloc_new_dim(&C, dim);
	int **A = (*matrix)[i1], **B = (*matrix)[i2];
	strassen_algorithm(A, B, &C, dim);
	return C;
}

// This function utilizes the usual multiplication algorithm and returns
// a pointer to the result.
int **usual_multiply(int ****matrix, int **m, int **n, int i1, int i2)
{
	int **mat;
	// Allocating memory
	mat = malloc((*m)[i1] * sizeof(int *));
	error_alloc2(mat);
	for (int i = 0; i < (*m)[i1]; ++i) {
		mat[i] = malloc((*n)[i2] * sizeof(int));
		error_alloc1(mat[i]);
	}

	for (int i = 0; i < (*m)[i1]; i++)
		for (int j = 0; j < (*n)[i2]; j++) {
			// Initializing the value of each element to 0
			mat[i][j] = 0;

			// Calculating each element
			for (int k = 0; k < (*m)[i2]; k++)
				mat[i][j] = (mat[i][j] + ((*matrix)[i1][i][k] *
					(*matrix)[i2][k][j]) + MOD * MOD) % MOD;
		}
	return mat;
}

// This function multiplies two specified matrixes together
// and adds the result to the vector.
void multiply(int ****matrix, int **m, int **n, int *count_matrix, char c)
{
	// Declaring their indexes and reading them
	int i1, i2;
	scanf("%d %d", &i1, &i2);

	// Checking them
	if (check_index(i1, *count_matrix) && check_index(i2, *count_matrix)) {
		// Checking if multiplication if possible
		if ((*n)[i1] == (*m)[i2]) {
			// Decalring a pointer to the new matrix
			int **mat;

			if (c == 'M') {
				// Utilising the usual algorithm
				mat = usual_multiply(matrix, m, n, i1, i2);
			} else {
				// Utilising the Strassen algorithm
				mat = strassen_multiply(matrix, m, i1, i2);
			}

			// Allocating more memory for the vector
			alloc_vector(matrix, m, n, *count_matrix);

			// Adding the new matrix and its dimensions
			(*m)[*count_matrix] = (*m)[i1];
			(*n)[*count_matrix] = (*n)[i2];
			(*matrix)[*count_matrix] = mat;

			// Adding one to the number of matrixes in the vector
			++(*count_matrix);
		} else {
			// Printing an error message
			printf("Cannot perform matrix multiplication\n");
		}
	}
}

// This function swaps two values.
void swap_two_values(int *a, int *b)
{
	int aux;
	aux = *b;
	*b = *a;
	*a = aux;
}

// This function swaps two matrixes in the vector,
// as well as their dimensions.
void swap_two_matrixes(int ****matrix, int **m, int **n, int i1, int i2)
{
	int **aux;
	aux = (*matrix)[i1];
	(*matrix)[i1] = (*matrix)[i2];
	(*matrix)[i2] = aux;

	swap_two_values(&(*m)[i1], &(*m)[i2]);
	swap_two_values(&(*n)[i1], &(*n)[i2]);
}

// This function sorts the matrixes based on the sum of all their elements.
void sort_vector(int ****matrix, int **m, int **n, int count_matrix)
{
	// Declaring and allocating memory for a temporary array that stores the
	// sum of each matrix
	int *S = malloc(count_matrix * sizeof(int));
	error_alloc1(S);

	// Initializing the value of each element to 0
	for (int i = 0; i < count_matrix; ++i)
		S[i] = 0;

	// Calucalting each sum
	for (int i = 0; i < count_matrix; ++i)
		for (int j = 0; j < (*m)[i]; ++j)
			for (int k = 0; k < (*n)[i]; ++k)
				S[i] = (S[i] + (*matrix)[i][j][k] + MOD) % MOD;

	// Sorting, using a selection sort algorithm
	for (int i = 0; i < count_matrix - 1; ++i)
		for (int j = i; j < count_matrix; ++j)
			if (S[i] > S[j]) {
				// Swaping the matrixes if necessary
				swap_two_matrixes(matrix, m, n, i, j);
				// Making sure to swap the sums as well
				swap_two_values(&S[i], &S[j]);
			}

	// Freeing the temporary array
	free(S);
}

// This function transposes a specified matrix.
void transpose(int ****matrix, int **m, int **n, int count_matrix)
{
	// Declaring its index and reading it
	int index;
	scanf("%d", &index);

	// Checking it
	if (check_index(index, count_matrix)) {
		// Declaring a pointer to the new matrix
		int **mat;

		// Allocating memory for it
		mat = malloc((*n)[index] * sizeof(int *));
		error_alloc2(mat);
		for (int i = 0; i < (*n)[index]; ++i) {
			mat[i] = malloc((*m)[index] * sizeof(int));
			error_alloc1(mat[i]);
		}

		// Creating the transpose
		for (int i = 0; i < (*n)[index]; ++i)
			for (int j = 0; j < (*m)[index]; ++j)
				mat[i][j] = (*matrix)[index][j][i];

		// Freeing the memeory allocated for the original matrix
		for (int i = 0; i < (*m)[index]; ++i)
			free((*matrix)[index][i]);
		free((*matrix)[index]);

		// Adding the transpose in the original's place
		(*matrix)[index] = mat;

		// Swaping the number of rows and columns
		swap_two_values(&(*m)[index], &(*n)[index]);
	}
}

// This function frees the memory allocated for a specified matrix.
// It also changes the arrays that hold the dimensions accordingly.
void free_matrix(int ****matrix, int **m, int **n, int *count_matrix)
{
	// Declaring its index and reading it
	int index;
	scanf("%d", &index);

	// Checking it
	if (check_index(index, *count_matrix)) {
		// Putting the specified matrix last in the vector
		for (int i = index + 1; i < *count_matrix; ++i) {
			swap_two_matrixes(matrix, m, n, index, i);
			index = i;
		}

		// Freeing the memory allocated for the last matrix
		// which is now the specified one
		for (int i = 0; i < (*m)[*count_matrix - 1]; ++i)
			free((*matrix)[*count_matrix - 1][i]);
		free((*matrix)[*count_matrix - 1]);

		// Subtracting one from the number of matrixes in the vector
		--(*count_matrix);

		// Allocating less memory to the vector and to the arrays
		alloc_vector(matrix, m, n, *count_matrix);
	}
}

// This function frees all memory allocated to the vector and to the arrays
// which store the dimensions of each matrix.
void free_vector(int ***matrix, int *m, int *n, int count_matrix)
{
	for (int i = 0; i < count_matrix; ++i) {
		for (int j = 0; j < m[i]; ++j)
			free(matrix[i][j]);
		free(matrix[i]);
	}
	free(matrix);
	free(m);
	free(n);
}

#endif
