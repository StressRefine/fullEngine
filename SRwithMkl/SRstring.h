/*
Copyright (c) 2020 Richard King

The stressRefine analysis executable "SRwithMkl" is free software:
you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SRwithMkl is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING.txt,
also available at <https://www.gnu.org/licenses/>
*/

//////////////////////////////////////////////////////////////////////
//
// SRstring.h: interface for the SRstring class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRSTRING_INCLUDED)
#define SRSTRING_INCLUDED

#include <string.h>
#include <string>
#include <vector>

using namespace std;


class SRstring
{
public:
	SRstring(SRstring& s2);
	SRstring(const char* s);
	SRstring();
	~SRstring();
	void Clear();
	bool isCsv();
	bool isAllBlank();
	bool caseInsensitiveCompare(const char* s2, int n = 0);
	void realStringCopy(char* dest, const char* src, int len);
	const char* getStr();
	const char* LastChar(const char c, bool after = false);
	void Copy(SRstring& s2);
	char GetChar(int i);
	char operator [] (int i);
	const char* FirstChar(char c);
	void operator = (SRstring& s2);
	void operator = (const char* s);
	void operator += (const char* s);
	void operator += (SRstring& s);
	bool operator == (SRstring& s2);
	bool Compare(SRstring& s2, int n = 0);
	bool CompareUseLength(const char* s2);
	bool CompareUseLength(SRstring& s2);
	bool operator == (const char* s2);
	bool operator != (const char* s2);
	void Left(char c, SRstring &s2, bool last = true);
	void Left(int n, SRstring &s2);
	void Right(char c, SRstring &s2);
	void Copy(const char* s, int n = 0);
	void Cat(const char* s);
	void Cat(SRstring& s2);
	bool Compare(const char* s2, int n = 0);
	bool isCommentOrBlank(bool skipContinuation = false);
	bool isBlank();
	bool isBdfComment(bool &isMat, SRstring& matname);
	const char* Token();
	const char* BdfToken(bool skipOnly = true);
	bool TokRead(int& i);
	bool TokRead(double& r, bool checkForComment = false);
	bool BdfRead(int& i);
	bool BdfRead(double& r);
	double RealRead();
	int IntRead();
	void setTokSep(const char sep);
	bool continueCheck();
	void TrimWhiteSpace();
	bool strIsBlank(const char* str);

	int FirstCharLocation(const char c);
	void bdfCheckLargeField();
	int getBdfWidth();
	void truncate(int n);
	int getLength();

	string str;
	int tokNum;
	char tokSep;
	int bdfPointer;
	int bdfWidth;
	bool fresh;
	vector <string> strSubs;
};

#endif //if !defined(SRSTRING_INCLUDED)
