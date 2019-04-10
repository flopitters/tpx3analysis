#include <iostream>
#include "hello.h"

using namespace std;


void hello_world ()
{
	cout << "Hello World!" << endl;
}

void print_message (char& str)
{
	cout << str << endl;
}

void duplicate (int& a, int& b, int& c)
{
	a *= 2;
	b *= 2;
	c *= 2;
}

void print_config (char& str)
{
	cout << str << endl;
}