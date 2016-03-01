#include <stdio.h>
#include <stdlib.h>

/**
 * Read a line of up to \size\ bytes from \fp\ to the buffer
 * @return 1 if done with file, 0 otherwise
 */ 
int readLine(FILE* fp, int size, char* buffer) {
	if (fp) {
		int pos = 0, c = EOL;
		while (pos < size-1) {
			c = fgetc(fp);
			if (c != EOL && c != '\n') {
				buffer[pos++] = (char) c;
			} else {
				break;
			}
		}
		buffer[pos] = (char) 0;
		if (c != EOL) return 0; // Not done
	}
	
	return 1; // Done
}

void trimEnd(int size, char* buffer) {
	// Find end of null-terminated string
	int pos = 0;
	while (buffer[pos] != 0 && pos < size-1) {
		pos += 1;
	}
	
	if (pos > 0 && buffer[pos] == 0) {
		while (buffer[pos-1] == ' ') {
			pos -= 1;
		}
		buffer[pos] = (char) 0;
	}
}

int main() {
	FILE* fp = fopen("equations.csv", "r");
	
	const int size = 64;
	char *buffer = (char*) malloc(size);
	
	int done = 0;
	do {
		done = readLine(fp, size, buffer);
		
	} while (done != 1);
	
	fclose(fp);
	free(buffer);
	
	return 0;
}