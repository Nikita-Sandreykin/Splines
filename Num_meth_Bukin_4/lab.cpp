#include "chartdir.h"
#include <vector>
#include <iostream>
#define pi atan(1)*4
#include <cmath>
using namespace std;

double f(double x)
{
	return 8/(pow((12+pi*x),3)) - pi*x;
}


class cubic_spline
{
	double tabX[4];
	double tabY[4];
	double a[3]; double b[3]; double c[3]; double d[3];
	double cAdd[4]; 
	double h[3]; double l[3];
	double q[2]; double v[2];
	public:
	void setTabPoints(double *tabX, double *tabY)
	{
		for (int i = 0; i < 4; i++)
		{
			this->tabX[i] = tabX[i];
			this->tabY[i] = tabY[i];
		}
	}
	void AC()
	{
		for (int i = 0; i < 3; i++)
		{
			*(h+i) = *(tabX + i + 1) - *(tabX+i);
			*(l+i) = (*(tabY + i + 1) - *(tabY+i)) / *(h+i);
		}
		cAdd[0] = 0;
	}
	void RC()
	{
		*q = -*(h+1) / (2 * (*h + *(h + 1)));
		*v = 3 * (*(l+1) - *l) / (2 * (*h + *(h+1)));
		for (int i = 2; i < 3; i++)
		{
			q[i - 1] = -h[i] / (2 * h[i - 1] + 2 * h[i] + h[i - 1] * q[i - 2]);
			v[i - 1] = (3 * l[i] - 3 * l[i - 1] - h[i - 1] * v[i - 2]) / (2 * h[i - 1] + 2 * h[i] + h[i - 1] * q[i - 2]);
		}
	
	}
	void SC()
	{
		cAdd[3] = 0;
		for (int i = 2; i > 0; i--)
		{
			cAdd[i] = q[i - 1] * cAdd[i + 1] + v[i - 1];
		}
		for (int i = 0; i < 3; i++)
		{
			a[i] = tabY[i + 1];
		}
		for (int i = 0; i < 3; i++)
		{
			c[i] = cAdd[i + 1];
		}
		for (int i = 0; i < 3; i++)
		{
			b[i] = l[i] + (2 * cAdd[i + 1] * h[i] + h[i] * cAdd[i]) / 3;
			d[i] = (cAdd[i + 1] - cAdd[i]) / (3 * h[i]);
		}
	}
	double getPoint(double x)
	{
		int k = 0;
		for (int i = 0; i < 3; i++)
		{
			if (x >= tabX[i] && x <= tabX[i + 1])
			{
				k = i;
			}
		}

		return a[k] + b[k] * (x - tabX[k + 1]) + c[k] * (x - tabX[k + 1]) * (x - tabX[k + 1]) + d[k] * (x - tabX[k + 1]) * (x - tabX[k + 1]) * (x - tabX[k + 1]);
	}
};
class parabolic_spline
{
private:
	double tabX[4];
	double tabY[4];
	double a[3]; 
	double b[3];
	double c[3];
public:
	void setTabPoints(double* tabX, double* tabY)
	{
		for (int i = 0; i < 4; i++)
		{
			this->tabX[i] = tabX[i];
			this->tabY[i] = tabY[i];
		}
	}
	void SC()
	{
		*c = 0; *(c + 2) = 0; 
		*b = (*tabY - *(tabY + 1)) / (*tabX - *(tabX + 1)); *(b + 2) = (*(tabY+2) - *(tabY+3)) / (*(tabX+2) - *(tabX+3));
		*(c+1) = (*b - *(b+2)) / (*(tabX+1) - *(tabX+2));
		*(b+1) = (*(tabY+1) - *(tabY+2) - *(c+1) * (*(tabX+1) * *(tabX + 1) - *(tabX+2) * *(tabX + 2))) / (*(tabX+1) - *(tabX+2));
		a[0] = tabY[1] - b[0] * tabX[1] - c[0] * tabX[1] * tabX[1]; a[1] = tabY[2] - b[1] * tabX[2] - c[1] * tabX[2] * tabX[2]; a[2] = tabY[3] - b[2] * tabX[3] - c[2] * tabX[3] * tabX[3];
	}
	double getPoint(double x)
	{
		int k = 0;
		for (int i = 0; i < 3; i++)
		{
			if (x >= tabX[i] && x <= tabX[i + 1])
			{
				k = i;
			}
		}
		return a[k] + b[k] * (x)+c[k] * (x) * (x);
	}
};
class linear_spline
{
	private:
	double tabX[4];
	double tabY[4];
	double a[3]; 
	double b[3];
public:
	void setTabPoints(double* tabX, double* tabY)
	{
		for (int i = 0; i < 4; i++)
		{
			this->tabX[i] = tabX[i];
			this->tabY[i] = tabY[i];
		}
	}
	void SC()
	{
		for (int i = 0; i < 3; i++)
		{
			*(a+i) = (*(tabY+i) - *(tabY + i + 1)) / (*(tabX+i) - *(tabX + i + 1));
			*(b+i) = *(tabY+i) - *(a+i) * *(tabX+i);
		}
	}
	double getPoint(double x)
	{
		int k = 0;
		for (int i = 0; i < 3; i++)
		{
			if (x >= tabX[i] && x <= tabX[i + 1])
			{
				k = i;
			}
		}
		return a[k] * x + b[k];
	}
};
int main(int argc, char* argv[])
{
	double X[] = { 2 , 2.5833, 3.1666, 3.75  };
	double Y[] = { f(2), f(2.5833), f(3.1666), f(3.75) };
	linear_spline spline1;
	parabolic_spline spline2;
	cubic_spline spline3;
	spline1.setTabPoints(X, Y); spline1.SC();
	spline2.setTabPoints(X, Y); spline2.SC();
	spline3.setTabPoints(X, Y); spline3.AC(); spline3.RC(); spline3.SC();
	vector<double> vfX, vlinearX, vparabolicX, vcubicX;
	vector<double> vfY, vlinearY, vparabolicY, vcubicY;
	vector<double> vel, vep, vec;
	for (double x = 2; x < 2.00001; x += 0.0000001)
	{
		vfX.push_back(x);  vlinearX.push_back(x); vparabolicX.push_back(x); vcubicX.push_back(x);
		vfY.push_back(f(x));  vlinearY.push_back(spline1.getPoint(x)); vparabolicY.push_back(spline2.getPoint(x)); vcubicY.push_back(spline3.getPoint(x));
		vel.push_back(abs(spline1.getPoint(x)-f(x))); vep.push_back(abs(spline2.getPoint(x) - f(x))); vec.push_back(abs(spline3.getPoint(x) - f(x)));
	}
	int n = vlinearX.size();
	cout << n << endl;
	cout << vfX.size();
	double *fX = new double[n];
	double *linearX = new double[n];
	double *parabolicX = new double[n];
	double *cubicX = new double[n];
	double* fY = new double[n];
	double *linearY = new double[n];
	double *parabolicY = new double[n];
	double *cubicY = new double[n];
	double* el = new double[n];
	double* ep = new double[n];
	double* ec = new double[n];
	for (int i = 0; i < n; i++)
	{
		fX[i] = vfX[i];
		fY[i] = vfY[i];
		linearX[i] = vlinearX[i];
		linearY[i] = vlinearY[i];
		parabolicX[i] = vparabolicX[i];
		parabolicY[i] = vparabolicY[i];
		cubicX[i] = vcubicX[i];
		cubicY[i] = vcubicY[i];
		el[i] = vel[i]; ep[i] = vep[i]; ec[i] = vec[i];
	}
	// Create a XYChart object of size 450 x 450 pixels
	XYChart* e = new XYChart(550, 550);
	XYChart* c = new XYChart(550, 550);
	c->setPlotArea(55, 65, 350, 300, 0xffffff, -1, 0xc0c0c0, 0xc0c0c0, -1);
	e->setPlotArea(55, 65, 350, 300, 0xffffff, -1, 0xc0c0c0, 0xc0c0c0, -1);
	c->addLegend(50, 30, false, "timesbi.ttf", 12)->setBackground(Chart::Transparent);
	e->addLegend(50, 30, false, "timesbi.ttf", 12)->setBackground(Chart::Transparent);
	c->addTitle("Splines", "timesbi.ttf", 18);
	e->addTitle("Splines errors", "timesbi.ttf", 18);
	c->yAxis()->setTitle("Y", "arialbi.ttf", 12);
	e->yAxis()->setTitle("Y", "arialbi.ttf", 12);
	c->yAxis()->setWidth(3);
	e->yAxis()->setWidth(3);
	c->xAxis()->setTitle("X", "arialbi.ttf", 12);
	e->xAxis()->setTitle("X", "arialbi.ttf", 12);
	c->xAxis()->setWidth(3);
	e->xAxis()->setWidth(3);
	LineLayer* layer1 = c->addLineLayer(DoubleArray(linearY, n), 0xff3333, "Linear Spline");
	layer1->setXData(DoubleArray(linearX, n));
	LineLayer* layer0 = c->addLineLayer(DoubleArray(fY, n), 0xf1f90c, "f(x)");
	layer0->setXData(DoubleArray(fX, n));
	layer1->setLineWidth(1);
	LineLayer* layer2 = c->addLineLayer(DoubleArray(parabolicY, n), 0x0678FC, "Parabolic Spline");
	layer2->setXData(DoubleArray(parabolicX, n));
	layer2->setLineWidth(1);
	LineLayer* layer3 = c->addLineLayer(DoubleArray(cubicY, n), 0x06FC60, "Cubic Spline");
	layer3->setXData(DoubleArray(cubicX, n));
	LineLayer* layerE1 = e->addLineLayer(DoubleArray(el, n), 0xff3333, "Line error Spline");
	layerE1->setXData(DoubleArray(fX, n));
	LineLayer* layerE2 = e->addLineLayer(DoubleArray(ep, n), 0x0678FC, "Parabolic error Spline");
	layerE2->setXData(DoubleArray(fX, n));
	LineLayer* layerE3 = e->addLineLayer(DoubleArray(ec, n), 0x06FC60, "Cubic error Spline");
	layerE3->setXData(DoubleArray(fX, n));
	layer3->setLineWidth(1);
	c->makeChart("splines.png");
	e->makeChart("errors.png");
	//free up resources
	delete c; delete e; delete[] fX; delete[] fY; delete[] el; delete[] ep; delete[] ec;
	delete[] linearX; delete[] linearY; delete[] parabolicX; delete[] parabolicY; delete[] cubicX; delete[] cubicY;
	return 0;
}

