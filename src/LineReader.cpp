/*
 * LineReader.cpp
 *
 *  Created on: 23/10/2013
 *      Author: jorge
 */

#include "LineReader.h"

LineReader::LineReader(){
	_size1 = DEFAULT_NUM_INDS/32;
	_size2 = DEFAULT_NUM_INDS/32;
	_length1 = 0;
	_length2 = 0;
	_internalPos1 = 0;
	_internalPos2 = 0;
	_buffer = new uint8_t[_size1];
	_buffer = new uint8_t[SNP_MAX_NAME_SIZE];
	if (_buffer == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_bufferCase0 = new uint32_t[_size1];
	if (_bufferCase0 == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_bufferCase1 = new uint32_t[_size1];
	if (_bufferCase1 == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_bufferCase2 = new uint32_t[_size1];
	if (_bufferCase2 == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_bufferCtrl0 = new uint32_t[_size2];
	if (_bufferCtrl0 == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_bufferCtrl1 = new uint32_t[_size2];
	if (_bufferCtrl1 == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_bufferCtrl2 = new uint32_t[_size2];
	if (_bufferCtrl2 == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}

	_fileBufferSentinel = 0;
	_fileBufferLength = 0;
	_fileBufferR = new uint8_t[4096 + 8];
	if (_fileBufferR == NULL) {
		Utils::exit("Memory allocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	_fileBuffer = _fileBufferR + 8; /*make it aligned*/

	_searchArrSize = 32;
	_searchArr = new int[_searchArrSize];
}

LineReader::~LineReader(){
	if (_buffer) {
		delete[] _buffer;
	}

	if (_bufferCase0) {
		delete[] _bufferCase0;
	}

	if (_bufferCase2) {
		delete[] _bufferCase2;
	}

	if (_bufferCtrl0) {
		delete[] _bufferCtrl0;
	}

	if (_bufferCtrl2) {
		delete[] _bufferCtrl2;
	}

	if(_searchArr) {
		delete[] _searchArr;
	}
}

int LineReader::readTFAMLine(MyFilePt& pt, uint32_t line){

	int ch = '\n';

	// Skip the empty lines
	while (ch == '\n'){
		ch = myfgetc(pt);
	}

	if(ch == -1){
		return -1;
	}

	// Skip the family information
	while ((ch = myfgetc(pt)) != -1 && !isblank(ch)){
	}
	if (ch == -1) { // It means that there is no value
		Utils::exit("TFAM file: Incomplete individual information in line %d\n", line);
	}

	// Skip the name
	while ((ch = myfgetc(pt)) != -1 && !isblank(ch)){
	}
	if (ch == -1) { // It means that there is no value
		Utils::exit("TFAM file: Incomplete individual information in line %d\n", line);
	}

	// Following information: father (space) mother (space) sex (space) case
	for(int i=0; i <=6; i++){
		if((ch = myfgetc(pt)) == -1){
			Utils::exit("TFAM file: Incomplete individual information in line %d\n", line);
		}
	}

	if(ch < '1' || ch > '2'){
		Utils::exit("TFAM file: %d is a wrong identification for line %d\n", ch, line);
	}

	int chEnd = myfgetc(pt);
	if((chEnd != '\n') && (chEnd != '\r')) // Value to pass to the next line
		Utils::exit("TFAM file: Too many values for line %d\n", line);

	if(ch == '1'){
		return 1;
	}

	return 0;
}


bool LineReader::readTPEDLine(MyFilePt& pt, SNP* readSNP, uint32_t line, uint32_t numInds,
		bool* indsClass){
	int ch;

	// Skip the empty lines
	while ((ch = myfgetc(pt)) == '\n'){
	}

	if (ch == -1) { // It means that there is no value
		return false;
	}

	// Skip the chromosome information
	while ((ch != -1) && !isblank(ch)){
		ch = myfgetc(pt);
	}
	if (ch == -1) { // It means that there is no value
		Utils::exit("TPED file: Incomplete SNP information in line %d\n", line);
	}

	// Read the SNP name (only one line)
	_length1 = 0;
	while ((ch = myfgetc(pt)) != -1 && !isblank(ch)){
		if(_length1 == SNP_MAX_NAME_SIZE){
			Utils::exit("TPED file: SNP name longer than %d in line %d\n", SNP_MAX_NAME_SIZE, line);
		}
		_buffer[_length1++] = ch;
	}
	if (ch == -1) { // It means that there is no value
		Utils::exit("TPED file: Incomplete SNP information in line %d\n", line);
	}
	_buffer[_length1] = '\0';

	// Trim characters /[12]$ in case that the name is too long for the line of the file
	if (_length1 > 2 && _buffer[_length1 - 2] == '/'
			&& (_buffer[_length1 - 1] == '1' || _buffer[_length1 - 1] == '2')) {
		_length1 -= 2;
		_buffer[_length1] = '\0';
	}

	// Save the SNP name
	readSNP->setName(_buffer);

	// Skip the marker id
	while ((ch = myfgetc(pt)) != -1 && !isblank(ch)){
	}
	if (ch == -1) { // It means that there is no value
		Utils::exit("TPED file: Incomplete SNP information in line %d\n", line);
	}

	// Skip the genetic distance
	while ((ch = myfgetc(pt)) != -1 && !isblank(ch)){
	}
	if (ch == -1) { // It means that there is no value
		Utils::exit("TPED file: Incomplete SNP information in line %d\n", line);
	}

	_length1 = 0;
	_length2 = 0;
	_internalPos1 = 0;
	_internalPos2 = 0;

	// To represent the base information of this SNP
	// The first base that appears is considered as major
	int readBases = 1;
	uint8_t phen;

	// Look which is the minor and major allele
	// Use an additional pointer
	if((ch = myfgetc(pt)) == -1)
		Utils::exit("TPED file: Not enough values for the SNP %d\n", line);
	int majorAllele = ch;
	int minorAllele = 0;
	_searchArr[0] = ch;

	while(!minorAllele && (readBases<(2*numInds))){ // Look for the other allele
		// Skip blank spaces
		ch = myfgetc(pt);
		while ((ch == ' ') || (ch == '\t')){
			ch = myfgetc(pt);
		}
		if(ch == -1)
			Utils::exit("TPED file: Not second allele SNP %d\n", line);

		if(ch != majorAllele){
			minorAllele = ch;
		}
		if(readBases >= _searchArrSize){
			_resizeSearchArr(_searchArrSize);
		}
		_searchArr[readBases] = ch;
		readBases++;
	}

	if(readBases%2){ // Take another base to complete the pair
		// Skip blank spaces
		ch = myfgetc(pt);
		while ((ch == ' ') || (ch == '\t')){
			ch = myfgetc(pt);
		}
		if(ch == -1)
			Utils::exit("TPED file: Not second allele SNP %d\n", line);

		if(readBases >= _searchArrSize){
			_resizeSearchArr(_searchArrSize);
		}
		_searchArr[readBases] = ch;
		readBases++;
	}

	if(readBases < (2*numInds)){
		if(majorAllele > minorAllele){
			majorAllele = minorAllele;
		}
	}

	int readInds;
	readBases/=2;

	// Take the values already read from the file
	for(readInds=0; readInds<readBases; readInds++){
		phen = 0;
		// Take the first allele
		if(_searchArr[2*readInds] != majorAllele){
			phen++;
		}

		if(_searchArr[2*readInds+1] != majorAllele){
			phen++;
		}

		if(!indsClass[readInds]){ // It is a case
			_insertCase(phen);
		}
		else{ // It is a control
			_insertCtrl(phen);
		}
	}

	// Read until the end
	while(readInds < numInds){
		// Skip blank spaces
		ch = myfgetc(pt);
		while ((ch == ' ') || (ch == '\t')){
			ch = myfgetc(pt);
		}
		if(ch == -1)
			Utils::exit("TPED file: Not enough values for the SNP %d\n", line);

		phen = 0;
		// Take the first allele
		if(ch != majorAllele){
			phen++;
		}

		// Skip blank spaces
		ch = myfgetc(pt);
		while ((ch == ' ') || (ch == '\t')){
			ch = myfgetc(pt);
		}
		if(ch == -1)
			Utils::exit("TPED file: Not enough values for the SNP %d\n", line);

		// Take the second allele
		if(ch != majorAllele){
			phen++;
		}

		if(!indsClass[readInds]){ // It is a case
			_insertCase(phen);
		}
		else{ // It is a control
			_insertCtrl(phen);
		}

		readInds++;
	}

	ch = myfgetc(pt);
	if ((ch != '\n') && (ch != '\r'))
		Utils::exit("TPED file: Too many values for the SNP %d (ch == %d)\n", line, ch);

	// Insert the values obtained in the new SNP
	readSNP->setCaseValues(_bufferCase0, _bufferCase1, _bufferCase2, _length1 + (_internalPos1>0));
	readSNP->setCtrlValues(_bufferCtrl0, _bufferCtrl1, _bufferCtrl2, _length2 + (_internalPos2>0));

	return true;
}

void LineReader::_insertCtrl(uint8_t phen){
	// Each value is saved as 1 bit
	if (_length2 >= _size2) {
		_resizeBufferCtrl(_size2);
	}

	if(!_internalPos2){
		_bufferCtrl0[_length2] = 0;
		_bufferCtrl1[_length2] = 0;
		_bufferCtrl2[_length2] = 0;
	}
	else{
		_bufferCtrl0[_length2] = _bufferCtrl0[_length2] << 1;
		_bufferCtrl1[_length2] = _bufferCtrl1[_length2] << 1;
		_bufferCtrl2[_length2] = _bufferCtrl2[_length2] << 1;
	}

	if(!phen){
		_bufferCtrl0[_length2]++;
	}
	else if(phen == 2){
		_bufferCtrl2[_length2]++;
	}
	else{
		_bufferCtrl1[_length2]++;
	}

	_internalPos2++;
	if(_internalPos2 == 32){
		_internalPos2 = 0;
		_length2++;
	}
}

void LineReader::_insertCase(uint8_t phen){
	// Each value is saved as 1 bit
	if (_length1 >= _size1) {
		_resizeBufferCases(_size1);
	}

	if(!_internalPos1){
		_bufferCase0[_length1] = 0;
		_bufferCase1[_length1] = 0;
		_bufferCase2[_length1] = 0;
	}
	else{
		_bufferCase0[_length1] = _bufferCase0[_length1] << 1;
		_bufferCase1[_length1] = _bufferCase1[_length1] << 1;
		_bufferCase2[_length1] = _bufferCase2[_length1] << 1;
	}

	if(!phen){
		_bufferCase0[_length1]++;
	}
	else if(phen == 2){
		_bufferCase2[_length1]++;
	}
	else{
		_bufferCase1[_length1]++;
	}

	_internalPos1++;
	if(_internalPos1 == 32){
		_internalPos1 = 0;
		_length1++;
	}
}

void LineReader::_resizeBufferCases(size_t nsize) {
	if (nsize < _size1) {
		return;
	}

	// Allocate a new buffer
	_size1 = nsize * 2;
	uint32_t* nbuffer = new uint32_t[_size1];
	if (!nbuffer) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer, _bufferCase0, _length1*sizeof(uint32_t));

	// Release the old buffer
	delete[] _bufferCase0;
	_bufferCase0 = nbuffer;

	uint32_t* nbuffer1 = new uint32_t[_size1];
	if (!nbuffer1) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer1, _bufferCase1, _length1*sizeof(uint32_t));

	// Release the old buffer
	delete[] _bufferCase1;
	_bufferCase1 = nbuffer1;

	uint32_t* nbuffer2 = new uint32_t[_size1];
	if (!nbuffer2) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer2, _bufferCase2, _length1*sizeof(uint32_t));

	// Release the old buffer
	delete[] _bufferCase2;
	_bufferCase2 = nbuffer2;
}

void LineReader::_resizeBufferCtrl(size_t nsize) {
	if (nsize < _size2) {
		return;
	}

	// Allocate a new buffer
	_size2 = nsize * 2;
	uint32_t* nbuffer = new uint32_t[_size2];
	if (!nbuffer) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer, _bufferCtrl0, _length2*sizeof(uint32_t));

	// Release the old buffer
	delete[] _bufferCtrl0;
	_bufferCtrl0 = nbuffer;

	uint32_t* nbuffer1 = new uint32_t[_size2];
	if (!nbuffer1) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer1, _bufferCtrl1, _length2*sizeof(uint32_t));

	// Release the old buffer
	delete[] _bufferCtrl1;
	_bufferCtrl1 = nbuffer1;

	uint32_t* nbuffer2 = new uint32_t[_size2];
	if (!nbuffer2) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer2, _bufferCtrl2, _size2*sizeof(uint32_t));

	// Release the old buffer
	delete[] _bufferCtrl2;
	_bufferCtrl2 = nbuffer2;
}

void LineReader::printMI(MyFilePt& pt, uint32_t id1, uint32_t id2, uint32_t id3, float MIValue){
	int error = fprintf(pt, "%u %u %u %f\n", id1, id2, id3, MIValue);

	if(error < 0)
		Utils::exit("Output file: File write failed with error %d at line %d in file %s\n", error, __LINE__,
					__FILE__);
}

void LineReader::_resizeSearchArr(size_t nsize) {
	if (nsize < _searchArrSize) {
		return;
	}

	// Allocate a new buffer
	_searchArrSize = nsize * 2;
	int* nbuffer = new int[_searchArrSize];
	if (!nbuffer) {
		Utils::exit("Memory reallocation failed in file %s in line %d\n",
				__FUNCTION__, __LINE__);
	}
	// Copy the old data
	memcpy(nbuffer, _searchArr, nsize*sizeof(int));

	// Release the old buffer
	delete[] _searchArr;
	_searchArr = nbuffer;
}

