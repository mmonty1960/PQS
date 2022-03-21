/*
Particles Qualification Software
version for SolarPACES diffusion

Author: Marco Montecchi
        ENEA-Casaccia
        marco.montecchi@enea.it

This file is part of PQSexpo.

    PQSexpo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation version 3 of the License

    PQSexpo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2012-2018 Marco Montecchi
*/
#ifndef PQS_H
#define PQS_H
 
#include "ui_PQS.h"
#include <minpack.h>
#include <cminpack.h>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <gsl/gsl_math.h>
 
class PQS : public QWidget, private Ui::PQS_DLG
{
    Q_OBJECT
 
public:
    PQS(QWidget *parent= nullptr);
    QMap<QString, QLineEdit*> idToLineEdit;
    void setWin(const std::string& _winname);

private:
    Ui::PQS_DLG *ui;

 
public slots:
    void getFileTw();
    void getFileRw();
    void getFileRws();
    void getFileRwRef1();
    void getFileRref1();
    void getFileRwRef2();
    void getFileRref2();
    void getFileZLine();
    void getFileMref1();
    void setBaseWL();
    void Method1();
    void calcRpart();
    void setRange();
    void saveMRh();
    void closeEvent ( QCloseEvent * event );
    void plotter(int iGraph,int RD);
    void refPlot();
    void recalc();
};
 
 
#endif
