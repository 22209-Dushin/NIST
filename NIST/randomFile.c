#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef unsigned char	BitSequence;
typedef unsigned char BYTE;

void generateRandomFile(const char *filename, size_t numBytes) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file for writing.\n");
        return;
    }

    srand(time(NULL));

    for (size_t i = 0; i < numBytes; ++i) {
        unsigned char byte = rand() % 256;  
        fwrite(&byte, sizeof(byte), 1, file);
    }

    fclose(file);
    printf("Random file generated with %zu bytes.\n", numBytes);
}

BitSequence	*epsilon;
void
readsHexDigitsInBinaryFormat(FILE *fp, unsigned int N) {  //количество бит 
	int		i, done, bitsRead;
	BYTE	buffer[4];
	
	if ( (epsilon = (BitSequence *) calloc(N, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		return;
	}
		do {
			if ( fread(buffer, sizeof(unsigned char), 4, fp) != 4 ) {
				printf("READ ERROR:  Insufficient data in file.\n");
				free(epsilon);
				return;
			}
			done = convertToBit(buffer, 32, N, &bitsRead);
		} while ( !done );
	free(epsilon);
}

int
convertToBit(BYTE *x, int xBitLength, int bitsNeeded, int *bitsRead)
{
	int		i, j, count, bit;
	BYTE	mask;

	count = 0;
	for ( i=0; i<(xBitLength+7)/8; i++ ) {
		mask = 0x80;
		for ( j=0; j<8; j++ ) {
			if ( *(x+i) & mask ) {
				bit = 1;
			}
			else {
				bit = 0;
			}
			mask >>= 1;
			epsilon[*bitsRead] = bit;
			(*bitsRead)++;
			if ( *bitsRead == bitsNeeded )
				return 1;
			if ( ++count == xBitLength )
				return 0;
		}
	}
	
	return 0;
}

int main() {
    int n = 1000000 * 8;

    generateRandomFile("data.bin", (size_t) (n / 8));

    return 0;
}

//gcc dffft.c  cephes.c matrix.c mainf.c -o main2