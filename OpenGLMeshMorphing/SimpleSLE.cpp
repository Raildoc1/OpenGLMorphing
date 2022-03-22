//#include <corecrt_malloc.h>
//#include <vcruntime_string.h>
//
//float* solveSLE(float** A, float* b, int n) {
//	float** AA = (float**)calloc(n, sizeof(float*));
//	for (size_t i = 0; i < n; i++)
//	{
//		AA[i] = (float*)calloc(n + 1, sizeof(float));
//		memcpy(AA[i], A[i], n * sizeof(float));
//		AA[i][n + 1] = -b[i];
//	}
//
//	for (int i = 0; i < n - 1; i++) {
//		for (int h = i + 1; h < n; h++)
//		{
//			float t = AA[h][i] / AA[i][i];
//			for (int j = 0; j <= n; j++)
//			{
//				AA[h][j] = AA[h][j] - t * AA[i][j];
//			}
//		}
//	}
//
//	float* result = (float*)calloc(n, sizeof(float));
//
//	for (int i = n - 1; i >= 0; i--)
//	{
//		result[i] = AA[i][n];
//		for (int j = n - 1; j > i; j--)
//		{
//			result[i] = result[i] - AA[i][j] * result[j];
//		}
//		result[i] = result[i] / AA[i][i];
//	}
//
//	return result;
//}