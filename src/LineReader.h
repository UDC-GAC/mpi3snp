/*
 * LineReader.h
 *
 *  Created on: 23/10/2013
 *      Author: jorge
 */

#ifndef LINEREADER_H_
#define LINEREADER_H_

#include "MyFile.h"
#include "SNP.h"

class LineReader {
public:
	LineReader();
	virtual ~LineReader();

	// 0 if it is a case
	// 1 if it is a control
	// -1 if it is the end of the file
	int readTFAMLine(MyFilePt& pt, uint32_t line);

	// Returns false if there is no more information
	bool readTPEDLine(MyFilePt& pt, SNP* readSNP, uint32_t line, uint32_t numInds,
			bool* indsClass);

	// To print the information of the interaction
	void printMI(MyFilePt& pt, uint32_t id1, uint32_t id2, uint32_t id3, float MIValue);

private:
	/*buffered file operations*/
	inline int myfgetc(MyFilePt file) {
		/*check the end-of-file*/
		if (_fileBufferSentinel >= _fileBufferLength) {
			/*re-fill the buffer*/
			_fileBufferSentinel = 0;
			/*read file*/
			_fileBufferLength = myfread(_fileBuffer, 1, 4096, file);
			if (_fileBufferLength == 0) {
				/*reach the end of the file*/
				if (myfeof(file)) {
					return -1;
				} else {
					Utils::exit("File reading failed in function %s line %d\n",
							__FUNCTION__, __LINE__);
				}
			}
		}
		/*return the current character, and increase the sentinel position*/
		return _fileBuffer[_fileBufferSentinel++];
	}
	inline int myungetc(int ch, MyFilePt file) {
		if (_fileBufferSentinel >= 0) {
			_fileBuffer[--_fileBufferSentinel] = ch;
		} else {
			Utils::log("Two consecutive ungetc operations occurred\n");
			return -1; /*an error occurred, return end-of-file marker*/
		}
		return ch;
	}

	void _insertCase(uint8_t phen);
	void _insertCtrl(uint8_t phen);

	void _resizeBufferCases(size_t nsize);
	void _resizeBufferCtrl(size_t nsize);

	void _resizeSearchArr(size_t nsize);

	uint8_t* _fileBufferR;
	uint8_t* _fileBuffer;
	int _fileBufferLength;
	int _fileBufferSentinel;

	// Array to save the values that are explored before knowing the two bases of the SNP
	int* _searchArr;
	int _searchArrSize;

	uint8_t* _buffer;
	uint32_t* _bufferCase0;
	uint32_t* _bufferCase1;
	uint32_t* _bufferCase2;
	uint32_t* _bufferCtrl0;
	uint32_t* _bufferCtrl1;
	uint32_t* _bufferCtrl2;
	// First value for cases and second value for ctrl
	// Keeps the number of SNP stored
	uint _length1;
	uint _length2;
	uint _internalPos1;
	uint _internalPos2;
	// Keeps the size of the array that is allocated
	uint _size1;
	uint _size2;
};

#endif /* LINEREADER_H_ */
