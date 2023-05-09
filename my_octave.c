// Mitran Andrei-Gabriel, 313CA
#include <stdio.h>
#include <stdlib.h>
#include "my_octave_func.h"

int main(void)
{
	// Declaring the vector and the arrays that store the dimensions
	int ***matrix = malloc(0), *m = malloc(0), *n = malloc(0);

	// Declaring the number of matrixes in the vector
	int count_matrix = 0;

	// Declaring and reading the variable that stores the command
	char c;
	scanf("%c", &c);

	// Executing until the command Q is given, which frees all
	// dynamically allocated memory
	while (c != 'Q') {
		// Making sure that the character is a letter
		if (c >= 'A' && c <= 'Z')
			// Executing each command accordingly
			switch (c) {
			case 'L':
				add_to_vector(&matrix, &m, &n, &count_matrix);
				break;

			case 'D':
				show_dim(m, n, count_matrix);
				break;

			case 'P':
				show(matrix, m, n, count_matrix);
				break;

			case 'C':
				resize(&matrix, &m, &n, count_matrix);
				break;

			case 'M':
				multiply(&matrix, &m, &n, &count_matrix, c);
				break;

			case 'O':
				sort_vector(&matrix, &m, &n, count_matrix);
				break;

			case 'T':
				transpose(&matrix, &m, &n, count_matrix);
				break;

			case 'F':
				free_matrix(&matrix, &m, &n, &count_matrix);
				break;

			case 'S':
				multiply(&matrix, &m, &n, &count_matrix, c);
				break;

			// Printing an error message if the command is unrecognized
			default:
				printf("Unrecognized command\n");
			}
		// Reading the variable again
		scanf("%c", &c);
	}

	// Freeing all dynamically allocated memory as per the Q command,
	// which was received earlier
	free_vector(matrix, m, n, count_matrix);

	return 0;
}
