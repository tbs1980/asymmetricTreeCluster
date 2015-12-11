#ifndef STEP_REGRESSION_HPP
#define STEP_REGRESSION_HPP

#include <assert.h>
#include <Eigen/Core>

class StepRegression
{
public:
	typedef Eigen::VectorXd VectorXd;
	typedef VectorXd::Scalar Scalar;
	typedef VectorXd::Index Index;

	StepRegression()
	{
		m_ready = false;
	}

	StepRegression(VectorXd const & xvals, VectorXd const & yvals)
	{
		setup(xvals,yvals);
	}

	~StepRegression()
	{

	}

	void setup(VectorXd const & xvals, VectorXd const & yvals)
	{
		// x,y data must be of the same size
		assert(xvals.rows()==yvals.rows());
		// There must be at least two data points
		assert(xvals.rows()>1);

		// Store x,y dats points
		m_x = xvals;
		m_y = yvals;
		m_N = xvals.rows();

		// Check data points are sorted in x
		for(Index i=1;i<m_N;i++) {assert(m_x(i) > m_x(i-1));}

		m_ready = true;
	}

	inline Scalar operator()(Scalar const& x) const
	{
		// Check x value is in bounds
		assert( x >= m_x(0) && x <= m_x(m_N-1));
		// There muse be at least two data points
		assert(m_N>1);
		// Check that the input data has been provided
		assert(m_ready == true);

		// Intial bounding solutions
		Index klo,khi,k;
		klo = 0;
		khi = m_N-1;

		// Binary search for lower and upper bounding x values
		while (khi-klo > 1)
		{
			k=(khi+klo) >> 1;
			if(m_x(k) > x) {khi=k;}
			else           {klo=k;}
		}

		// Assert that the solutions are within the range of input x values
		assert(khi>=0 && klo>=0);
		assert(khi<m_N && klo<m_N);

		// Step regression: use the low-bounding x value to predict y
		return m_y(klo);
	}

private:

	VectorXd m_x;
	VectorXd m_y;
	Index m_N;
	bool m_ready;
};

#endif // STEP_REGRESSION_HPP
