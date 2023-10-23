#include <math.h>

double TaskOne(int x)
{
	double ans = 0;
	int end = std::max(20, 20 * abs(x));
	for (int k = 1; k <= end; k++)
	{
		for (int j = 1; j <= end; j++)
		{
			ans += ((x * x * (k - j)) / (x * x + pow(k, 3) + pow(j, 3))) * sin(k * x) * cos(j * x);
		}
	}
	return ans;
}