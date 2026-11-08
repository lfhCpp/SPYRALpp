#ifndef matrix_op_h
#define matrix_op_h


namespace matrix_op
{

using size_t = std::size_t;
#ifdef USE_SINGLE_PRECISION
using Real_T = float;
#else
using Real_T = double;
#endif


/* `replace_row<OPT:N_ROW, M_COL>(row_idx, row, mat)`:
 * Replaces a row of an matrix by a vector.
 */
template <size_t N_ROW, size_t M_COL>
static void
replace_row(const size_t row_idx, const Real_T (&row)[M_COL], Real_T (&mat)[N_ROW][M_COL])
{
	for (size_t i = 0; i < M_COL; ++i) {
		mat[row_idx][i] = row[i];
	}
}


/*  `weighted_sum(u_weight, u, v_weight, v, OUT:u_weight * u + v_weight * v)`:
 * Sums two vectors after weighting.
 */
template <size_t DIM>
void
weighted_sum(const Real_T u_weight, const Real_T (&u)[DIM], const Real_T v_weight,
             const Real_T (&v)[DIM], Real_T (&sum)[DIM])
{
	for (size_t i = 0; i < DIM; ++i) {
		sum[i] = u_weight * u[i] + v_weight * v[i];
	}
}

/*  `weighted_sum(A_weight, A, B_weight, B, OUT:A_weight * A + B_weight * B)`:
 * Sums two matrices after weighting.
 */
template <size_t N_DIM, size_t M_DIM>
void
weighted_sum(const Real_T A_weight, const Real_T (&A)[N_DIM][M_DIM], const Real_T B_weight,
             const Real_T (&B)[N_DIM][M_DIM], Real_T (&sum)[N_DIM][M_DIM])
{
	for (size_t i = 0; i < N_DIM; ++i) {
		for (size_t j = 0; j < M_DIM; ++j) {
			sum[i][j] = A_weight * A[i][j] + B_weight * B[i][j];
		}
	}
}



} // namespace matrix_op



#endif