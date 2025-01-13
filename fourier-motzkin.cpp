#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//first column is b, then returns vector x such that b <= A*x
//cf. https://people.math.carleton.ca/~kcheung/math/notes/MATH5801/02/2_1_fourier_motzkin.html#FM
vec fourier_motzkin (const mat& A, const double tol=1e-6, const double the_best_number = 1) {
	if (A.n_elem == 0) {
		return vec(A.n_cols-1, fill::value(the_best_number));
	}

	vec a = A.col(A.n_cols-1); // eliminate last column
	a.clean(tol);
	uvec lower_bounds = find(a > 0);
	uvec upper_bounds = find(a < 0);
	uvec carries      = find(a == 0);
	if (A.n_cols == 1) {
		if (lower_bounds.n_elem > 0) throw logic_error("system is infeasible");
		else return vec();
	}

	a(carries).ones(); //leave carries alone

	mat smaller_A (carries.n_elem+lower_bounds.n_elem*upper_bounds.n_elem, A.n_cols);
	smaller_A.head_rows(carries.n_elem) = A.rows(carries);
	uword i = carries.n_elem;
	for (uword l : lower_bounds) for (uword u : upper_bounds)
		smaller_A.row(i++) = A.row(l)/a(l) - A.row(u)/a(u);
	smaller_A.shed_col(smaller_A.n_cols-1); //remove last column which is zero now

	vec smaller_x = fourier_motzkin(smaller_A);

	vec b = A.col(0);
	if (A.n_cols > 2)
		b -= A.cols(1, A.n_cols-2) * smaller_x;
	b /= a;

	const double upper = upper_bounds.n_elem>0 ? b.rows(upper_bounds).min() : nan("upper");
	const double lower = lower_bounds.n_elem>0 ? b.rows(lower_bounds).max() : nan("lower");
	const double less_than_zero = carries.n_elem>0 ? b.rows(carries).max()  : 0;
	if ( tol+upper < lower || tol < less_than_zero)
		throw logic_error("system is infeasible");

	double x;
	if (isnan(upper) && isnan(lower)) x = the_best_number; //unbounded in both directions
	else if (isnan(upper) && !isnan(lower)) x = lower;
	else if (!isnan(upper) && isnan(lower)) x = upper;
	else x = upper;
	smaller_x.insert_rows(smaller_x.n_rows, vec(1, fill::value(x)));

	return smaller_x;
}

int main() {
	//mat A = {{0,1,1},{2,2,1},{1,-1,1},{-1,-1,2}};
	mat A = {{2,1,1,-2},{0,-1,-3,1},{1,0,1,1}};
	//mat A = {{1,1,1,2},{2,-1,1,1},{1,1,-1,1},{0,0,-1,-3}};
	cout << A << endl;

	try {
		vec x = fourier_motzkin(A);
		cout << "Result: " << endl << x << endl;
	} catch(const logic_error& e) {
		cout << e.what() << endl;
	}
	return 0;
}
