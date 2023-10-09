#include <iostream>
#include "Random64.h"

int main(void)
{
	Crandom rand64(1);
	for(int i=0;i<100;i++)
	{
		std::cout << rand64.r()-1 << std::endl;
	}
	
	return 0;
}