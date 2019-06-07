#include "colormapgenerator.h"

#include <limits>
#include <cmath>

ColorMapGenerator::ColorMapGenerator(bool useZeroAsMin)
{
    if (useZeroAsMin)
        this->setFixedMinValue(0.0f);
    else
        this->useDynamicMinValue();

    this->useDynamicMaxValue();
    this->setColorMapToHeat();
}

QVector<QColor> ColorMapGenerator::createColorMap(QVector<BSplineSurface *> bSplines, BSplineTriangulation *triangulation, SubdivisionAbstractMetric *metric)
{
    const int numSamplePoints = triangulation->getNumberOfVertices();

    QVector<double> metricValues(numSamplePoints);

    #pragma omp parallel for
    for (int i = 0; i < numSamplePoints; i++) {
        const int surfId = triangulation->getSurfaceIdOfVertex(i);
        const QVector2D param = triangulation->getParameterValueForVertex(i);
        metricValues[i] = metric->evaluatePoint(bSplines.at(surfId), param.x(), param.y());
    }

    //min/max is not thread safe -> own for loop
    double min = metricValues[0];
    double max = metricValues[0];
    for (int i = 1; i < numSamplePoints; i++) {
        const double value = metricValues[i];
        if (value < min)
            min = value;
        if (value > max)
            max = value;
    }

    if (!m_useFixedMaximum)
        m_maxValue = max;

    if (!m_useFixedMinimum)
        m_minValue = min;

    const double valueRange = m_maxValue - m_minValue;

    QVector<QColor> colorMap(numSamplePoints);
    #pragma omp parallel for
    for (int i = 0; i < numSamplePoints; i++) {
        double relativeValue = (metricValues[i] - m_minValue) / valueRange;

        if (relativeValue <= 0.0f) {
            colorMap[i] = m_colorMapColors.first();
        } else if (relativeValue >= 1.0f) {
            colorMap[i] = m_colorMapColors.last();
        } else {
            int interval = 0;
            while (relativeValue > m_colorMapPercentiles[interval + 1]) //value < 1 = colormapPercentiles.last() (due to above if-check)
                interval++;

            const double lower = m_colorMapPercentiles[interval];
            const double upper = m_colorMapPercentiles[interval + 1];
            const double intervalWidth = upper - lower;
            const double lambda = (relativeValue - lower)/intervalWidth;    //(1 - lambda) * lower + (lambda * upper) = relativeValue

            const QVector3D col0(m_colorMapColors[interval].red(), m_colorMapColors[interval].green(), m_colorMapColors[interval].blue());
            const QVector3D col1(m_colorMapColors[interval + 1].red(), m_colorMapColors[interval + 1].green(), m_colorMapColors[interval + 1].blue());
            const QVector3D col = (1 - lambda) * col0 + lambda * col1;

            colorMap[i] = QColor(floor(col.x() + 0.5), floor(col.y() + 0.5), floor(col.z() + 0.5));
        }
    }

    return colorMap;
}

double ColorMapGenerator::getMinMetricValue()
{
    return m_minValue;
}

double ColorMapGenerator::getMaxMetricValue()
{
    return m_maxValue;
}

QString ColorMapGenerator::generateMatlabPlotCommand(const int colorMapSamples, const QString cBarFilename)
{
    QString result("cm = [");
    for (int i = 0; i < colorMapSamples; i++) {

        const double relativeValue = (double) i /(colorMapSamples - 1);
        int interval = 0;
        while (relativeValue > m_colorMapPercentiles[interval + 1]) //value <= 1 = colormapPercentiles.last()
            interval++;

        const double lower = m_colorMapPercentiles[interval];
        const double upper = m_colorMapPercentiles[interval + 1];
        const double intervalWidth = upper - lower;
        const double lambda = (relativeValue - lower)/intervalWidth;    //(1 - lambda) * lower + (lambda * upper) = relativeValue

        const QVector3D col0(m_colorMapColors[interval].red(), m_colorMapColors[interval].green(), m_colorMapColors[interval].blue());
        const QVector3D col1(m_colorMapColors[interval + 1].red(), m_colorMapColors[interval + 1].green(), m_colorMapColors[interval + 1].blue());
        const QVector3D col = (1 - lambda) * col0 + lambda * col1;

        result.append(QString("%1 %2 %3").arg(col.x()/255).arg(col.y()/255).arg(col.z()/255));
        if (i < colorMapSamples - 1)
            result.append("; ");
    }
    result.append("];\n");
    result.append(QString("plotColorBar(%1, %2, cm, '%3');").arg(m_minValue).arg(m_maxValue).arg(cBarFilename));

    return result;
}

void ColorMapGenerator::useDynamicMaxValue()
{
    m_useFixedMaximum = false;
}

void ColorMapGenerator::useDynamicMinValue()
{
    m_useFixedMinimum = false;
}

void ColorMapGenerator::setFixedMaxValue(double maxValue)
{
    m_maxValue = maxValue;
    m_useFixedMaximum = true;
}

void ColorMapGenerator::setFixedMinValue(double minValue)
{
    m_minValue = minValue;
    m_useFixedMinimum = true;
}

void ColorMapGenerator::setColorMapToHeat(bool includeWhite)
{
    m_colorMapColors.resize(4);
    m_colorMapPercentiles.resize(4);
    m_colorMapColors[0] = QColor(0,   0,   0  );
    m_colorMapColors[1] = QColor(255, 0,   0  );
    m_colorMapColors[2] = QColor(255, 255, 0  );
    if (includeWhite)
        m_colorMapColors[3] = QColor(255, 255, 255);
    else
        m_colorMapColors[3] = QColor(255, 255, 200);
    m_colorMapPercentiles[0] = 0.0f;
    m_colorMapPercentiles[1] = 0.2f;
    m_colorMapPercentiles[2] = 0.6f;
    m_colorMapPercentiles[3] = 1.0f;
}

void ColorMapGenerator::setCustomColorMap(const QVector<QColor> &colors, const QVector<double> &indices)
{
    if (m_colorMapColors.size() != m_colorMapPercentiles.size())
        return;
    m_colorMapColors = colors;
    m_colorMapPercentiles = indices;
}
