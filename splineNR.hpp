#ifndef SPLINE_NR_HPP
#define SPLINE_NR_HPP

#include <assert.h>
#include <Eigen/Core>

class gen1DSpline
{
public:
	typedef Eigen::VectorXd VectorXd;
	typedef VectorXd::Scalar Scalar;
	typedef VectorXd::Index Index;

	gen1DSpline()
	// :m_x(0),m_y(0),m_N(0),m_ready(false)
	:m_ready(false)
	{

	}

	gen1DSpline(VectorXd const & xvals, VectorXd const & yvals)
	{
		assert(xvals.rows()==yvals.rows());
		assert(xvals.rows()>2);

		interpolate(xvals,yvals);
	}

	~gen1DSpline()
	{

	}

	void interpolate(VectorXd const & xvals, VectorXd const & yvals)
	{
		assert(xvals.rows()==yvals.rows());
		assert(xvals.rows()>2);

		m_x = xvals;
		m_y = yvals;
		m_N = xvals.rows();
		m_y2 = VectorXd(m_N);

		//check if ordered
		for(Index i=1;i<m_N;i++)
		{
			assert(m_x(i) > m_x(i-1));
		}

		VectorXd u(m_N);
		Scalar qn,un;

		m_y2(0)=0;
		u(0)=0;

		for (Index i=1;i<m_N-1;i++)
		{
			Scalar sig = ( m_x(i)-m_x(i-1) )/( m_x(i+1)-m_x(i-1) );
			Scalar p = sig*m_y2(i-1) + 2.;
			m_y2(i) = (sig-1.)/p;
			u(i) = ( m_y(i+1)-m_y(i) )/( m_x(i+1)-m_x(i) )
				- ( m_y(i)-m_y(i-1) )/( m_x(i)-m_x(i-1) );
			u(i) = ( 6.*u(i)/( m_x(i+1)-m_x(i-1) )-sig*u(i-1))/p;
		}

		qn = 0.;
		un = 0.;

		m_y2(m_N-1) = ( un-qn*u(m_N-2) )/( qn*m_y2(m_N-2)+1. );



		for(Index k=m_N-1;k>0;--k)//This is a problem?
		{
			m_y2(k-1) = m_y2(k-1)*m_y2(k) + u(k-1);
		}

		m_ready = true;


	}

	inline Scalar operator()(Scalar x) const
	{
		assert( x >= m_x(0) && x <= m_x(m_N-1));
		//assert(m_N>2);
		assert(m_ready == true);

		Index klo,khi,k;
		Scalar h,b,a,y;
		klo=0;
		khi=m_N-1;
		while (khi-klo > 1)
		{
			k=(khi+klo) >> 1;
			if(m_x(k) > x)
			{
				khi=k;
			}
			else
			{
				klo=k;
			}
		}

		assert(khi<m_N && klo<m_N);

		h = m_x(khi)-m_x(klo);

		assert(std::abs(h) > 0.);

		a = ( m_x(khi)-x )/h;
		b = ( x-m_x(klo) )/h;

		y = a*m_y(klo) + b*m_y(khi)+( (a*a*a-a)*m_y2(klo) + (b*b*b-b)*m_y2(khi) )*(h*h)/6.;

		return y;
	}

	inline Index nPoints() const
	{
		return m_N;
	}

	inline Scalar xMax(void) const
	{
		return m_x(m_N-1);
	}

	inline Scalar xMin(void) const
	{
		return m_x(0);
	}

	void writeToTextFile(std::string const& fileName) const
	{
		Index precision=10;
		std::ofstream file;
		file.open(fileName.c_str(),std::ios::trunc);

		if(file.is_open())
		{
			file<<std::scientific;
			file<<std::setprecision(precision);

			for(Index i=0;i<m_N;++i)
			{
				file<<m_x(i)<<","<<m_y(i)<<std::endl;
			}

			file.close();
		}
		else
		{
			std::string msg = "Error opening file:\t"	+ fileName ;
			throw std::runtime_error(msg);
		}

	}

private:

	VectorXd m_x;
	VectorXd m_y;
	VectorXd m_y2;
	Index m_N;
	bool m_ready;

};

#endif //SPLINE_NR_HPP
