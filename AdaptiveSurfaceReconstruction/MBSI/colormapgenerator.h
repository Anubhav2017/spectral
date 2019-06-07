#ifndef COLORMAPGENERATOR_H
#define COLORMAPGENERATOR_H

#include "BSpline/bsplinesurface.h"
#include "BSpline/bsplinetriangulation.h"
#include "subdivisionabstractmetric.h"

#include <QString>
#include <QColor>

class ColorMapGenerator
{
public:
    ColorMapGenerator(bool useZeroAsMin);

    QVector<QColor> createColorMap(QVector<BSplineSurface *> bSplines, BSplineTriangulation *triangulation, SubdivisionAbstractMetric *metric);
    double getMinMetricValue();
    double getMaxMetricValue();
    QString generateMatlabPlotCommand(const int colorMapSamples, const QString cBarFilename);

    void useDynamicMaxValue();
    void useDynamicMinValue();
    void setFixedMaxValue(double maxValue);
    void setFixedMinValue(double minValue);
    void setColorMapToHeat(bool includeWhite = true);
    void setCustomColorMap(const QVector<QColor> &colors, const QVector<double> &indices);
private:
    QVector<QColor> m_colorMapColors;
    QVector<double> m_colorMapPercentiles;

    bool m_useFixedMaximum;
    bool m_useFixedMinimum;
    double m_minValue;
    double m_maxValue;
};

#endif // COLORMAPGENERATOR_H
