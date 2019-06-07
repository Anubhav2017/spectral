/********************************************************************************
** Form generated from reading UI file 'ImageInterface.ui'
**
** Created by: Qt User Interface Compiler version 5.9.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMAGEINTERFACE_H
#define UI_IMAGEINTERFACE_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_ImageInterface
{
public:
    QVBoxLayout *vboxLayout;
    QGridLayout *gridLayout;
    QSpinBox *imgHeight;
    QSpinBox *imgWidth;
    QLabel *label_2;
    QSpacerItem *spacerItem;
    QLabel *label_3;
    QCheckBox *ratioCheckBox;
    QHBoxLayout *hboxLayout;
    QLabel *label;
    QDoubleSpinBox *oversampling;
    QSpacerItem *spacerItem1;
    QComboBox *color_comboBox;
    QCheckBox *expandFrustum;
    QSpacerItem *spacerItem2;
    QHBoxLayout *hboxLayout1;
    QSpacerItem *spacerItem3;
    QPushButton *okButton;
    QPushButton *cancelButton;

    void setupUi(QDialog *ImageInterface)
    {
        if (ImageInterface->objectName().isEmpty())
            ImageInterface->setObjectName(QStringLiteral("ImageInterface"));
        ImageInterface->resize(484, 447);
        vboxLayout = new QVBoxLayout(ImageInterface);
        vboxLayout->setSpacing(6);
        vboxLayout->setObjectName(QStringLiteral("vboxLayout"));
        vboxLayout->setContentsMargins(9, 9, 9, 9);
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        imgHeight = new QSpinBox(ImageInterface);
        imgHeight->setObjectName(QStringLiteral("imgHeight"));
        imgHeight->setMinimum(1);
        imgHeight->setMaximum(32000);

        gridLayout->addWidget(imgHeight, 0, 4, 1, 1);

        imgWidth = new QSpinBox(ImageInterface);
        imgWidth->setObjectName(QStringLiteral("imgWidth"));
        imgWidth->setMinimum(1);
        imgWidth->setMaximum(32000);

        gridLayout->addWidget(imgWidth, 0, 1, 1, 1);

        label_2 = new QLabel(ImageInterface);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 0, 0, 1, 1);

        spacerItem = new QSpacerItem(20, 22, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(spacerItem, 0, 2, 1, 1);

        label_3 = new QLabel(ImageInterface);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout->addWidget(label_3, 0, 3, 1, 1);

        ratioCheckBox = new QCheckBox(ImageInterface);
        ratioCheckBox->setObjectName(QStringLiteral("ratioCheckBox"));
        ratioCheckBox->setChecked(true);

        gridLayout->addWidget(ratioCheckBox, 1, 4, 1, 1);


        vboxLayout->addLayout(gridLayout);

        hboxLayout = new QHBoxLayout();
        hboxLayout->setSpacing(6);
        hboxLayout->setObjectName(QStringLiteral("hboxLayout"));
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(ImageInterface);
        label->setObjectName(QStringLiteral("label"));

        hboxLayout->addWidget(label);

        oversampling = new QDoubleSpinBox(ImageInterface);
        oversampling->setObjectName(QStringLiteral("oversampling"));
        oversampling->setDecimals(1);
        oversampling->setMinimum(0.1);
        oversampling->setMaximum(64);
        oversampling->setSingleStep(1);
        oversampling->setValue(1);

        hboxLayout->addWidget(oversampling);

        spacerItem1 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacerItem1);


        vboxLayout->addLayout(hboxLayout);

        color_comboBox = new QComboBox(ImageInterface);
        color_comboBox->setObjectName(QStringLiteral("color_comboBox"));

        vboxLayout->addWidget(color_comboBox);

        expandFrustum = new QCheckBox(ImageInterface);
        expandFrustum->setObjectName(QStringLiteral("expandFrustum"));

        vboxLayout->addWidget(expandFrustum);

        spacerItem2 = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Expanding);

        vboxLayout->addItem(spacerItem2);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setSpacing(6);
        hboxLayout1->setObjectName(QStringLiteral("hboxLayout1"));
        hboxLayout1->setContentsMargins(0, 0, 0, 0);
        spacerItem3 = new QSpacerItem(131, 31, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout1->addItem(spacerItem3);

        okButton = new QPushButton(ImageInterface);
        okButton->setObjectName(QStringLiteral("okButton"));

        hboxLayout1->addWidget(okButton);

        cancelButton = new QPushButton(ImageInterface);
        cancelButton->setObjectName(QStringLiteral("cancelButton"));

        hboxLayout1->addWidget(cancelButton);


        vboxLayout->addLayout(hboxLayout1);


        retranslateUi(ImageInterface);
        QObject::connect(okButton, SIGNAL(clicked()), ImageInterface, SLOT(accept()));
        QObject::connect(cancelButton, SIGNAL(clicked()), ImageInterface, SLOT(reject()));

        QMetaObject::connectSlotsByName(ImageInterface);
    } // setupUi

    void retranslateUi(QDialog *ImageInterface)
    {
        ImageInterface->setWindowTitle(QApplication::translate("ImageInterface", "Image Settings", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        imgHeight->setToolTip(QApplication::translate("ImageInterface", "Height of the image (in pixels)", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        imgHeight->setSuffix(QApplication::translate("ImageInterface", " px", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        imgWidth->setToolTip(QApplication::translate("ImageInterface", "Width of the image (in pixels)", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        imgWidth->setSuffix(QApplication::translate("ImageInterface", " px", Q_NULLPTR));
        label_2->setText(QApplication::translate("ImageInterface", "Width", Q_NULLPTR));
        label_3->setText(QApplication::translate("ImageInterface", "Height", Q_NULLPTR));
        ratioCheckBox->setText(QApplication::translate("ImageInterface", "Keep Ratio", Q_NULLPTR));
        label->setText(QApplication::translate("ImageInterface", "Oversampling", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        oversampling->setToolTip(QApplication::translate("ImageInterface", "Antialiases image (when larger then 1.0)", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        oversampling->setPrefix(QApplication::translate("ImageInterface", "x ", Q_NULLPTR));
        color_comboBox->clear();
        color_comboBox->insertItems(0, QStringList()
         << QApplication::translate("ImageInterface", "Use current background color", Q_NULLPTR)
         << QApplication::translate("ImageInterface", "Use transparent background", Q_NULLPTR)
         << QApplication::translate("ImageInterface", "Choose background color", Q_NULLPTR)
        );
#ifndef QT_NO_TOOLTIP
        expandFrustum->setToolTip(QApplication::translate("ImageInterface", "When image aspect ratio differs from viewer's one, expand frustum as needed. Fits inside current frustum otherwise.", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        expandFrustum->setText(QApplication::translate("ImageInterface", "Expand Frustum if Needed", Q_NULLPTR));
        okButton->setText(QApplication::translate("ImageInterface", "OK", Q_NULLPTR));
        cancelButton->setText(QApplication::translate("ImageInterface", "Cancel", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ImageInterface: public Ui_ImageInterface {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMAGEINTERFACE_H
