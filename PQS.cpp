/*
Paricle Qualification software
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

    Copyright 2012-2019 Marco Montecchi
*/
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <qtextstream.h>
#include "PQS.h"
#include <minpack.h>
#include <cminpack.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <qwt_plot_grid.h>
#include <qwt_legend.h>
#include <qwt_symbol.h>
#include <qwt_interval.h>
#include <qwt_plot_intervalcurve.h>
#include <qwt_interval_symbol.h>
#include <gsl/gsl_math.h>
#include <tgmath.h>

using namespace std;

//global variables *******************************
QString pathroot;
QString pathroot2;
QString fileTwnk,fileRwnk,fileZLnk,fileRefnk;
QString fileTw,fileRw,fileRws,fileRwRef1,fileRef1,fileZL,fileRwRef2,fileRef2,fileMref1;
int Nw;
int iUVN1IR2=1;
double chi2fin,A,B,C,xy[2][100],RhExp[5000][5],dW;
double Pig=3.14159265359;

int iBGR;//index Rdiff column data

QwtPlot *G_Rpart;
QColor myColor[7]={Qt::black,Qt::blue,Qt::cyan,Qt::green,Qt::magenta,Qt::red,Qt::gray};
int IXW=-1;//plot window not yet opened

double TEMAs,TEMAp,TEMAu,REMAs,REMAp,REMAu,nOl,kOl;
double MRh[3000][30];
/*
  MRh[3000][101] main data storage matrix
        [i][0] = wavelength (nm)
        [i][1] = weight (absolute)
        [i][2] =
        [i][3] =
        [i][4] =
        [i][5] =
        [i][6] = n_Window
        [i][7] = k_Window
        [i][8] = <Tw>
        [i][9] = <Rw>
        [i][10] =
        [i][11] = Tw
        [i][12] = Rw
        [i][13] = Rws
        [i][14] = RwRef1
        [i][15] = Rref1
        [i][16] = RwRef2
        [i][17] = Rref2
        [i][18] = ZLine
        [i][19] = MeasRef1
        [i][20] = Rpart_Meth0
        [i][21] = Rpart_Meth1
        [i][22] = Rpart_Meth2
        [i][23] = Rpart_Meth3
        [i][24] = Rpart_Meth4
        [i][25] = Rpart_Meth5
        [i][26] =
        [i][27] =
        [i][28] =
        [i][29] =
*/


// funzioni invocate
void parabolicFit(int n);
void spada(int column, QString fileRh, double division,int IUVN1IR2);
void RsRpEMA(int iWL, double theta);
void RsRp(complex<double> N1, complex<double> N2, complex<double> theta1,
          double& Ts,double& Tp,double& Rs,double& Rp, complex<double>& theta2);
void nkSUB(double lambda, double T, double R);


PQS::PQS(QWidget *parent){
    setupUi(this); // this sets up GUI
    Q_UNUSED( parent )

    // signals/slots mechanism in action
    connect(pushButton_Tw,  SIGNAL(clicked()), this, SLOT(getFileTw()));
    connect(pushButton_Rw,  SIGNAL(clicked()), this, SLOT(getFileRw()));
    connect(pushButton_Rws, SIGNAL(clicked()), this, SLOT(getFileRws()));
    connect(pushButton_RwRef1, SIGNAL(clicked()), this, SLOT(getFileRwRef1()));
    connect(pushButton_Rref1, SIGNAL(clicked()), this, SLOT(getFileRref1()));
    connect(pushButton_ZLine,  SIGNAL(clicked()), this, SLOT(getFileZLine()));
    connect(pushButton_Mref1, SIGNAL(clicked()), this, SLOT(getFileMref1()));
    connect(pushButton_RwRef2, SIGNAL(clicked()), this, SLOT(getFileRwRef2()));
    connect(pushButton_Rref2, SIGNAL(clicked()), this, SLOT(getFileRref2()));
    connect(comboBox_standard,SIGNAL(currentIndexChanged(int)), this, SLOT(setRange()));
    connect(pushButton_calcRpart,SIGNAL(clicked()), this, SLOT(calcRpart()));
    connect(dSB_Ymin,SIGNAL(valueChanged(double)), this, SLOT(calcRpart()));
    connect(dSB_Ymax,SIGNAL(valueChanged(double)), this, SLOT(calcRpart()));
    connect(dSB_Lamb,SIGNAL(valueChanged(double)), this, SLOT(recalc()));
    connect(pushButton_refresh,SIGNAL(clicked()), this, SLOT(refPlot()));

#ifdef __unix__
#define IS_POSIX 1
#else
#define IS_POSIX 0
#endif

    if(IS_POSIX == 1) {
        //Linux path initialization
        const QByteArray value = qgetenv("USER");
        QString uName=QString::fromLocal8Bit(value);
        cout << "current user = " << uName.toStdString() <<endl;
        pathroot="/home/"+uName+"/Workspace/Particles";
    }
    else{
        //windows path inizialization
        pathroot=getenv("PWD");
        pathroot=pathroot+"workspace/Particles";
    }
    pathroot2=pathroot;
    setRange();//initialization
}

void PQS::closeEvent ( QCloseEvent *  ){
    qApp->quit();
}

void PQS::setRange(){
    int iStd = comboBox_standard -> currentIndex();
    comboBox_range -> clear();
    if(iStd<5){
        comboBox_range -> addItem("solar 320-2500 nm");
        comboBox_range -> addItem("UVA 315-400 nm");
        if(iStd==0 || iStd==1){
            comboBox_range -> addItem("UVB 280-315 nm");
            comboBox_range -> addItem("UV  280-400 nm");
        }
        iUVN1IR2=1;
    }
    else{
        comboBox_range -> addItem("IR 2.5-15.38 um");
        iUVN1IR2=2;
    }
}

void PQS::getFileTw(){
    fileTw = QFileDialog::getOpenFileName(
                this,
                "Choose the window-Transmittance file",
                pathroot2);
    if(fileTw.isEmpty())
        return;
    lineEdit_Tw -> setText(fileTw.section('/',-1));
}

void PQS::getFileRw(){
    fileRw = QFileDialog::getOpenFileName(
                this,
                "Choose the window-Reflectance file",
                pathroot2);
    if(fileRw.isEmpty())
        return;
    lineEdit_Rw -> setText(fileRw.section('/',-1));
}

void PQS::getFileRws(){
    fileRws = QFileDialog::getOpenFileName(
                this,
                "Choose the window-particle file",
                pathroot);
    if(fileRws.isEmpty())
        return;
    lineEdit_Rws -> setText(fileRws.section('/',-1));
    pathroot2=fileRws.section('/',1,-1);
}

void PQS::getFileRwRef1(){
    fileRwRef1 = QFileDialog::getOpenFileName(
                this,
                "Choose the window-Reference#1 file",
                pathroot2);
    if(fileRwRef1.isEmpty())
        return;
    lineEdit_RwRef1 -> setText(fileRwRef1.section('/',-1));
}

void PQS::getFileRref1(){
    fileRef1 = QFileDialog::getOpenFileName(
                this,
                "Choose the R_reference#1 file",
                pathroot+"/Standard");
    if(fileRef1.isEmpty())
        return;
    lineEdit_Rref1 -> setText(fileRef1.section('/',-1));
}

void PQS::getFileMref1(){
    fileMref1 = QFileDialog::getOpenFileName(
                this,
                "Choose the Meas_reference#1 file",
                pathroot2);
    if(fileMref1.isEmpty())
        return;
    lineEdit_Mref1 -> setText(fileMref1.section('/',-1));
}

void PQS::getFileRwRef2(){
    fileRwRef2 = QFileDialog::getOpenFileName(
                this,
                "Choose the window-Reference#2 file",
                pathroot2);
    if(fileRwRef2.isEmpty())
        return;
    lineEdit_RwRef2 -> setText(fileRwRef2.section('/',-1));
}

void PQS::getFileRref2(){
    fileRef2 = QFileDialog::getOpenFileName(
                this,
                "Choose the R_reference#2 file",
                pathroot+"/Standard");
    if(fileRef2.isEmpty())
        return;
    lineEdit_Rref2 -> setText(fileRef2.section('/',-1));
}

void PQS::getFileZLine(){
    fileZL = QFileDialog::getOpenFileName(
                this,
                "Choose the ZeroLine file",
                pathroot2);
    if(fileZL.isEmpty())
        return;
    lineEdit_ZLine -> setText(fileZL.section('/',-1));
}

void PQS::setBaseWL(){
    int iRg,iStd,iRowB=0,iRowE=0;
    QString fileW;

    //standard
    iRg=comboBox_range -> currentIndex();     //range
    iStd=comboBox_standard -> currentIndex(); //standard
    if(iStd==0)
        fileW=pathroot+"/Standard/IEC60904b3Step5nm.txt";
    else if(iStd==1)
        fileW=pathroot+"/Standard/ASTMG173SP.txt";
    else if(iStd==2)
        fileW=pathroot+"/Standard/ISO9050.txt";
    else if(iStd==3)
        fileW=pathroot+"/Standard/ISO9845b1.txt";
    else if(iStd==4)
        fileW=pathroot+"/Standard/E891.txt";
    else if(iStd==5)
        fileW=pathroot+"/Standard/IRnoWeight.txt";
    if(iStd<=1){
        if(iRg==0){
            iRowB=9;
            iRowE=445;
        }
        else if(iRg==1){
            iRowB=8;
            iRowE=25;
        }
        else if(iRg==2){
            iRowB=1;
            iRowE=8;
        }
        else if(iRg==3){
            iRowB=1;
            iRowE=25;
        }
    }
    else if(iStd<5){
        if(iRg==0){
            iRowB=5;
            iRowE=441;
        }
        else if(iRg==1){
            iRowB=4;
            iRowE=21;
        }
    }
    else{
        iRowB=1;
        iRowE=336;
    }
    QFile fW(fileW);
    if (!fW.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream stream0 (&fW);
    Nw=0;
    int iRow=0;
    double dum1,dum2;
    do {
        iRow++;
        if(iRowB<=iRow && iRow<=iRowE){
            stream0 >> MRh[Nw][0] >> MRh[Nw][1];
            //printf("%d	%f	%f \n", Nw,MRh[Nw][0],MRh[Nw][1]);
            Nw++;
        }
        else {
            stream0 >> dum1 >> dum2;
        }
    } while (!stream0.atEnd());
    printf("Loaded %s \n Nw=%d \n",(fileW.toStdString()).c_str(),Nw);
    fW.close();
}


void PQS::Method1(){//Window characterizatio for Method#1
    setBaseWL();

    int iStd=comboBox_standard -> currentIndex(); //standard
    dW=doubleSpinBox_dW->value();
    if(iStd!=5)
        dW=dW*1.e+06;//mm->nm
    else
        dW=dW*1.e+03;//mm->um
    double lambda,T,R,Rref1,ZL,ZLg;
    double expLamb=dSB_Lamb->value();
    //computing nk-window from Tw & Rw
    for(int iW=0;iW<Nw;iW++){
        lambda=MRh[iW][0];
        Rref1=MRh[iW][15];
        ZL=MRh[iW][18];
        ZLg=0.;
        if(iUVN1IR2==2){
            ZLg=ZL;
            ZL=0.;
        }
        T=MRh[iW][11]-ZLg;
        R=(MRh[iW][12]-ZLg-T*T*ZL)*Rref1;
        nkSUB(lambda, T, R);
        MRh[iW][6]=nOl;
        MRh[iW][7]=kOl;
        printf("wl=%f T=%f R=%f A=%f n=%f k=%e\n",lambda,T,R,1.-T-R,nOl,kOl);

        //computing <Tw> & <Rw>
        double sTw=0.,sRw=.0,theta;
        for(int j=0;j<90;j++){
            theta=double(j)/180*Pig;
            RsRpEMA(iW,theta);
            sTw=sTw+TEMAu*pow(cos(theta),expLamb)*sin(theta);
            sRw=sRw+REMAu*pow(cos(theta),expLamb)*sin(theta);
        }
        MRh[iW][8]=sTw*2.*Pig/180.;
        MRh[iW][9]=sRw*2.*Pig/180.;
        printf("WL=%f <Tw>=%f <Rw>=%f Tw-<Tw>=%f Rw-<Rw>=%f\n",
               MRh[iW][0],MRh[iW][8],MRh[iW][9],T-MRh[iW][8],R-MRh[iW][9]);
    }
}


void PQS::calcRpart(){
    setBaseWL();
    int iMeth=comboBox_Meth->currentIndex();
    if(iMeth==1){
        Method1();
    }
    printf("->calRpart iMeth=%d\n",iMeth);
    Qt::CheckState state1;
    state1 = checkBox -> checkState();
    double division=1.;
    if( state1 == Qt::Checked ) division=100.;
    //Tw
    spada(11,fileTw,division,iUVN1IR2);
    //printf("Loaded %s\n",(fileTw.toStdString()).c_str());
    //Rw
    spada(12,fileRw,division,iUVN1IR2);
    //printf("Loaded %s\n",(fileRw.toStdString()).c_str());
    //Rws
    spada(13,fileRws,division,iUVN1IR2);
    //printf("Loaded %s\n",(fileRws.toStdString()).c_str());
    //RwRef1
    spada(14,fileRwRef1,division,iUVN1IR2);
    //printf("Loaded %s\n",(fileRwRef1.toStdString()).c_str());
    //Rref1
    spada(15,fileRef1,1.,1);
    //printf("Loaded %s\n",(fileRef1.toStdString()).c_str());
    //RwRef2
    spada(16,fileRwRef2,division,iUVN1IR2);
    //printf("Loaded %s\n",(fileRwRef2.toStdString()).c_str());
    //Rref2
    spada(17,fileRef2,1.,1);
    //printf("Loaded %s\n",(fileRef2.toStdString()).c_str());
    //ZLine
    spada(18,fileZL,division,iUVN1IR2);
    //printf("Loaded %s\n",(fileZL.toStdString()).c_str());
    //MeasRef1
    if(iUVN1IR2==2)
        spada(19,fileMref1,division,iUVN1IR2);

    long double Tw,Rw,Rws,RwRef1,Rref1,RwRef2,Rref2,ZL,ZLg,TwM,RwM,Rpart,Norm,sumR=0.,sumW=0.;
    for(int iW=0;iW<Nw;iW++){
        Rref1=MRh[iW][15];
        Rref2=MRh[iW][17];
        ZL=MRh[iW][18];
        ZLg=0.;
        Norm=1.;
        if(iUVN1IR2==2){
            ZLg=ZL;
            ZL=0.;
            Norm=MRh[iW][19]-ZLg;
        }
        Tw=MRh[iW][11]-ZLg;
        Rw=(MRh[iW][12]-ZLg-Tw*Tw*ZL)*Rref1/Norm;
        Rws=(MRh[iW][13]-ZLg)*Rref1/Norm;
        RwRef1=(MRh[iW][14]-ZLg)*Rref1/Norm;
        RwRef2=(MRh[iW][16]-ZLg)*Rref2/Norm;
        if(iMeth==0){
            Rpart=(Rws-Rw)/(Tw*Tw-Rw*Rw+Rw*Rws);
            MRh[iW][20]=Rpart;
        }
        else if(iMeth==1){
            TwM=MRh[iW][8];
            RwM=MRh[iW][9];
            Rpart=(Rws-Rw)/(Tw*TwM+RwM*(Rws-Rw));
            MRh[iW][21]=Rpart;
        }
        else if(iMeth==2){
            long double x=(RwRef1-Rw)/Rref1*
                    (1.-Rref1/(RwRef2-RwRef1)*((RwRef2-Rw)/Rref2-(RwRef1-Rw)/Rref1));
            long double y=1./(RwRef2-RwRef1)*((RwRef2-Rw)/Rref2-(RwRef1-Rw)/Rref1);
            Rpart=(Rws-Rw)/(x+y*(Rws-Rw));
            MRh[iW][22]=Rpart;
        }
        else if(iMeth==3){
            TwM=(RwRef1-Rw)*(1-Rref1)/(Rref1*(Tw+Rw-RwRef1));
            Rpart=(Rws-Rw)/(Tw*TwM+(1.-TwM)*(Rws-Rw));
            MRh[iW][23]=Rpart;
            //printf("Meth_3a wl=%f <Tw>=%f Tw=%f\n",MRh[iW][0],double(TwM),double(Tw));
        }
        else if(iMeth==4){
            TwM=(1.-Rref1*Rw)*(RwRef1-Rw)/(Rref1*Tw);
            Rpart=(Rws-Rw)/(Tw*TwM+(1.-TwM)*(Rws-Rw));
            MRh[iW][24]=Rpart;
        }
        else{//5
            RwM=(RwRef1-Rw-Rref1*Tw*Tw)/(Rref1*(RwRef1-Rw));
            Rpart=(Rws-Rw)/(Tw*Tw+RwM*(Rws-Rw));
            MRh[iW][25]=Rpart;
        }
        if(1.-Tw-Rw < 0.1){
            sumW=sumW+MRh[iW][1];
            sumR=sumR+Rpart*MRh[iW][1];
        }
    }
    printf("<Rparticle>=%f\n",double(sumR/sumW));
    lineEdit_meanRp->setText(QString::number(double(sumR/sumW)));
    saveMRh();
    plotter(iMeth,0);
    plotter(10,0);
}

void PQS::recalc(){
    comboBox_Meth->setCurrentIndex(1);
    calcRpart();
}

void PQS::saveMRh(){
    QString fnSave=pathroot+"/MatrixRparticles.txt";
    QFile fS(fnSave);
    if (!fS.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out(&fS);
    out << "WL(nm)"<<"\t"<<"Rp_Meth#0"<<"\t"<<"Rp_Meth#1"<<"\t"<<"Rp_Meth#2"
        <<"\t"<<"Rp_Meth#3"<<"\t"<<"Rp_Meth#4"<<"Rp_Meth#5"<<"\n";
    for(int i=0; i<Nw; i++){
        out<<MRh[i][0]<<"\t"<<MRh[i][20]<<"\t"<<MRh[i][21]<<"\t"<<MRh[i][22]<<"\t"<<MRh[i][23]
                     <<"\t"<<MRh[i][24]<<"\t"<<MRh[i][25]<<"\n";
    }
    fS.close();

    QString fnSaveMRh=pathroot+"/MRh.txt";
    QFile fS2(fnSaveMRh);
    if (!fS2.open(QIODevice::WriteOnly | QIODevice::Text))
        return;
    QTextStream out2(&fS2);
    out2 << "WL(nm)"<<"\t"<<"Tw"<<"\t"<<"<Tw>"<<"\t"<<"Rw"<<"\t"<<"<Rw>"<<"\n";
    double Tw,TwM,Rw,RwM,ZL,ZLg,Rref1;
    for(int iW=0; iW<Nw; iW++){
        Rref1=MRh[iW][15];
        ZL=MRh[iW][18];
        ZLg=0.;
        if(iUVN1IR2==2){
            ZLg=ZL;
            ZL=0.;
        }
        Tw=MRh[iW][11]-ZLg;
        Rw=(MRh[iW][12]-ZLg-Tw*Tw*ZL)*Rref1;
        TwM=MRh[iW][8];
        RwM=MRh[iW][9];
        out2<<MRh[iW][0]<<"\t"<<Tw<<"\t"<<TwM<<"\t"<<Rw<<"\t"<<RwM<<"\n";
    }
    fS2.close();
}

void PQS::refPlot(){
    int iMeth=comboBox_Meth->currentIndex();
    plotter(iMeth,1);
}

void PQS::plotter(int iGraph,int iRD){
    // iGraph==0 -> plot Rpart Meth#0
    // iGraph==1 -> plot Rpart Meth#1
    // iGraph==2 -> plot Rpart Meth#2
    // iGraph==3 -> plot Rpart Meth#3a
    // iGraph==4 -> plot Rpart Meth#3b
    // iGraph==5 -> plot Rpart Meth#4
    // iGraph==10-> plot Aw
    // int iRD= 0->overwrite 1->redraw
    printf("->plotter iGraph=%d iRD=%d\n",iGraph,iRD);
    int iCol=0;
    double xP[Nw],yP[Nw];
    double WLmin=1.e+6;
    double WLmax=-1.;
    double Ymin=dSB_Ymin->value();
    double Ymax=dSB_Ymax->value();
    double Rref1,ZL,ZLg,T,R,Norm;
    for(int iW=0;iW<Nw;iW++){
        if(iGraph<10){
            xP[iW]=MRh[iW][0];
            yP[iW]=MRh[iW][20+iGraph];
            iCol=iGraph;
        }
        else if(iGraph==10){
            xP[iW]=MRh[iW][0];
            Rref1=MRh[iW][15];
            ZL=MRh[iW][18];
            ZLg=0.;
            Norm=1.;
            if(iUVN1IR2==2){
                ZLg=ZL;
                ZL=0.;
                Norm=MRh[iW][19]-ZLg;
            }
            T=MRh[iW][11]-ZLg;
            R=(MRh[iW][12]-ZLg-T*T*ZL)*Rref1/Norm;
            yP[iW]=1.-T-R;
            iCol=6;
        }
        WLmin=min(xP[iW],WLmin);
        WLmax=max(xP[iW],WLmax);
    }
    if(IXW<0){
        G_Rpart=new QwtPlot();
        IXW=1;
    }
    else{
        if(iRD==1){
            G_Rpart->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G_Rpart->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
    }
    QwtPlotCurve *dataPlot=new QwtPlotCurve("dataPlot");
    QwtPlotIntervalCurve *range_plot = new QwtPlotIntervalCurve("range");
    dataPlot->setSamples(xP,yP,Nw);
    dataPlot->setPen(QPen(myColor[iCol],2,Qt::SolidLine));
    dataPlot->attach(G_Rpart);
//    /* error bars */
//    QVector<QwtIntervalSample> range(Nw);
//    for(int i = 0; i < Nw; i++) {
//        range[i] = QwtIntervalSample(xP[i],yP[i]-ErrYp[i],yP[i]+ErrYp[i]);
//    }
//    QwtIntervalSymbol *errorbar = new QwtIntervalSymbol(QwtIntervalSymbol::Bar);
//    errorbar->setPen(QPen(myColor[iCol],1));
//    errorbar->setWidth(0);
//    range_plot->setSamples(range);
//    range_plot->setSymbol(errorbar);
//    range_plot->setStyle(QwtPlotIntervalCurve::NoCurve);
    range_plot->attach(G_Rpart);
    G_Rpart -> setAxisTitle(0,"Rparticles");
    if(iUVN1IR2==1)
        G_Rpart -> setAxisTitle(2,"wavelength (nm)");
    else
        G_Rpart -> setAxisTitle(2,"wavelength (um)");
    G_Rpart->setAxisScale(0,Ymin,Ymax,0);
    G_Rpart->setAxisScale(2,WLmin,WLmax,0);
    G_Rpart->setAutoReplot();
    G_Rpart->show();

    // Make the grid on
    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen(QPen(Qt::gray, 0.0, Qt::DotLine));
    grid->enableX(true);
    grid->enableXMin(true);
    grid->enableY(true);
    grid->enableYMin(true);
    grid->attach(G_Rpart);
    G_Rpart->replot();
    G_Rpart->show();
}


void parabolicFit(int n){
    long double matrix[3][4], ratio, a;
    long double sum_x=0.,sum_y=0.,sum_x2=0.,sum_x3=0.,sum_x4=0.,sum_xy=0.,sum_x2y=0.;
    int i, j , k;
    for(i = 0; i < n; i++){
        sum_x += xy[0][i];
        sum_y += xy[1][i];
        sum_x2 += pow(xy[0][i], 2);
        sum_x3 += pow(xy[0][i], 3);
        sum_x4 += pow(xy[0][i], 4);
        sum_xy += xy[0][i]*xy[1][i];
        sum_x2y += pow(xy[0][i], 2) * xy[1][i];
    }
    matrix[0][0] = n;
    matrix[0][1] = sum_x;
    matrix[0][2] = sum_x2;
    matrix[0][3] = sum_y;
    matrix[1][0] = sum_x;
    matrix[1][1] = sum_x2;
    matrix[1][2] = sum_x3;
    matrix[1][3] = sum_xy;
    matrix[2][0] = sum_x2;
    matrix[2][1] = sum_x3;
    matrix[2][2] = sum_x4;
    matrix[2][3] = sum_x2y;
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            if(i!=j){
                ratio = matrix[j][i]/matrix[i][i];
                for(k = 0; k < 4; k++){
                    matrix[j][k] -= ratio * matrix[i][k];
                }
            }
        }
    }
    for(i = 0; i < 3; i++){
        a = matrix[i][i];
        for(j = 0; j < 4; j++){
            matrix[i][j] /= a;
        }
    }
    C=matrix[0][3];
    B=matrix[1][3];
    A=matrix[2][3];
}


void spada(int ic, QString fileRh, double division, int IUVN1IR2){
    int i,Ndata,Nflag,Nflag0,n,istep,iAlarmMIN=0,iAlarmMAX=0,iPlot11=0;
    double Wa,Wb,WL,Yfin,DET,xPlot[10000],yPlot[10000];//err
    QString line,msg;
    QStringList list;
    if(fileRh.isEmpty())
        return;
    printf("spada ic=%d division=%f IUVN1IR2=%d file=",ic,division,IUVN1IR2);
    cout<<fileRh.toStdString()<<"\n";
    
    QFile fRh(fileRh);
    if (!fRh.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    QTextStream stream ( &fRh );
    Ndata=0;
    double val[5];
    int Narg=2;
    int Nskip=0;
    do{
        line = fRh.readLine();
        line=line.simplified();
        //cout<<line.toStdString();
        list=line.split(QRegExp("\\s+"));
        int nV=list.size();
        bool OK=false;
        if(nV==Narg){
            for(int index=0;index<nV;index++){
                val[index]=list.at(index).toDouble(&OK);
                if(!OK)
                    index=nV;
            }
            //cout<<"   ->>>>> "<<nV<<"\t"<<OK<<"\n";
            if(OK){
                for (int k=0;k<nV;k++){
                    RhExp[Ndata][k]=val[k];
                    //cout<<val[k]<<"\t";
                }
                //cout<<"\n";
                if(IUVN1IR2==2){
                    RhExp[Ndata][0]=1/RhExp[Ndata][0]*1.e4;//cm^-1 -> um
                    //printf("RhExp[%d][0]=%f RhExp[%d][1]=%f\n",Ndata,RhExp[Ndata][0],Ndata,RhExp[Ndata][1]);
                }
                Ndata++;
            }
            else
                Nskip++;
        }
        else{
            Nskip++;
            //cout<<"   ->>>>> "<<nV<<"\t"<<OK<<"\n";
        }
    }while(!stream.atEnd());

    printf("Ndati=%d Nskip_line=%d\n",Ndata,Nskip);
    fRh.close();
    if(RhExp[0][0] < RhExp[Ndata-1][0]){
        istep=1;
        Nflag0=0;
        if(RhExp[0][0] > MRh[0][0]+0.1){
            iAlarmMIN=1;
            printf("istep=%d RhExp[0][0]=%f MRh[0][0]=%f\n",
                   istep,RhExp[0][0],MRh[0][0]);
        }
        if(RhExp[Ndata-1][0] < MRh[Nw-1][0]-0.1){
            printf("istep=%d RhExp[Ndata-1][0]=%f MRh[Nw-1][0]=%f\n",
                   istep,RhExp[Ndata-1][0],MRh[Nw-1][0]);
            iAlarmMAX=1;
        }
    }
    else{
        istep=-1;
        Nflag0=Ndata-1;
        if(RhExp[Ndata-1][0] > MRh[0][0]+0.1){
            printf("istep=%d RhExp[Ndata-1][0]=%f MRh[0][0]=%f\n",
                   istep,RhExp[Ndata-1][0],MRh[0][0]);
            iAlarmMIN=1;
        }
        if(RhExp[0][0] < MRh[Nw-1][0]-0.1){
            printf("istep=%d RhExp[0][0]=%f MRh[Nw-1][0]=%f\n",
                   istep,RhExp[0][0],MRh[Nw-1][0]);
            iAlarmMAX=1;
        }
    }
    if(iAlarmMIN!=0){
        msg="WLmin > "+QString::number(MRh[0][0]);
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
        return;
    }
    if(iAlarmMAX!=0){
        msg="WLmax < "+QString::number(MRh[Nw-1][0]);
        QMessageBox msgBox;
        msgBox.setText(msg);
        msgBox.exec();
        return;
    }
    // resampling RhExp like weight
    int jmax=1;
    for(int j=1;j<=jmax;j++){
        Nflag=Nflag0;
        for(i=0;i<Nw;i++){
            WL=MRh[i][0];
            if(i>0)
                Wa=WL-(MRh[i][0]-MRh[i-1][0])/2.;
            else
                Wa=WL-(MRh[i+1][0]-MRh[i][0])/2.;
            if(i<Nw-1)
                Wb=WL+(MRh[i+1][0]-MRh[i][0])/2.;
            else
                Wb=WL+(MRh[i][0]-MRh[i-1][0])/2.;
            n=-1;
            do{
                if(RhExp[Nflag][0] >= Wa){
                    n++;
                    xy[0][n]=RhExp[Nflag][0];
                    xy[1][n]=RhExp[Nflag][j];
                }
                Nflag=Nflag+istep;
            } while(RhExp[Nflag][0] <= Wb && Nflag>=0 && Nflag<= Ndata-1);
            //      printf("Wa=%f Wb=%f n=%d \n",Wa,Wb,n);
            Nflag=Nflag-istep;
            Yfin=0.;
            if(n==-1){
                n=0;
                xy[0][n]=RhExp[Nflag][0];
                xy[1][n]=RhExp[Nflag][j];
                n=1;
                Nflag=Nflag+istep;
                xy[0][n]=RhExp[Nflag][0];
                xy[1][n]=RhExp[Nflag][j];
            }
            if(n==0){
                //	printf("Nflag=%d\n",Nflag);
                if(xy[0][n]<WL){
                    if(Nflag+istep >=0 && Nflag+istep <= Ndata-1){
                        //	    printf("caso 1\n");
                        Nflag=Nflag+istep;
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                    else{
                        //	    printf("caso 2\n");
                        n=0;
                        xy[0][n]=RhExp[Nflag-istep][0];
                        xy[1][n]=RhExp[Nflag-istep][j];
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                }
                else{
                    if(Nflag-istep >=0 && Nflag-istep <= Ndata-1){
                        //	    printf("caso 3\n");
                        n=0;
                        xy[0][n]=RhExp[Nflag-istep][0];
                        xy[1][n]=RhExp[Nflag-istep][j];
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                    else{
                        //	    printf("caso 4\n");
                        Nflag=Nflag+1;
                        n=1;
                        xy[0][n]=RhExp[Nflag][0];
                        xy[1][n]=RhExp[Nflag][j];
                    }
                }
            }
            if(n==1){
                Yfin=xy[1][0]+(xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0])*(WL-xy[0][0]);
            }
            if(n==2){
                DET=xy[0][0]*xy[0][0]*(xy[0][1]-xy[0][2])-xy[0][0]*(xy[0][1]*xy[0][1]-xy[0][2]*xy[0][2]);
                DET=DET+xy[0][1]*xy[0][1]*xy[0][2]-xy[0][2]*xy[0][2]*xy[0][1];
                A=xy[1][0]*(xy[0][1]-xy[0][2])-xy[0][0]*(xy[1][1]-xy[1][2]);
                A=(A+xy[1][1]*xy[0][2]-xy[1][2]*xy[0][1])/DET;
                B=xy[0][0]*xy[0][0]*(xy[1][1]-xy[1][2])-xy[1][0]*(xy[0][1]*xy[0][1]-xy[0][2]*xy[0][2]);
                B=(B+xy[0][1]*xy[0][1]*xy[1][2]-xy[0][2]*xy[0][2]*xy[1][1])/DET;
                C=xy[0][0]*xy[0][0]*(xy[0][1]*xy[1][2]-xy[0][2]*xy[1][1]);
                C=C-xy[0][0]*(xy[0][1]*xy[0][1]*xy[1][2]-xy[0][2]*xy[0][2]*xy[1][1]);
                C=(C+xy[1][0]*(xy[0][1]*xy[0][1]*xy[0][2]-xy[0][2]*xy[0][2]*xy[0][1]))/DET;
                Yfin=A*WL*WL+B*WL+C;
            }
            if(n>2){
                parabolicFit(n);
                Yfin=A*WL*WL+B*WL+C;
            }
//            for(int jj=0;jj<=n;jj++)
//                printf("%f %f\n",xy[0][jj],xy[1][jj]);
//            if(n>=2) printf("A=%f B=%f C=%f\n",A,B,C);
//            printf("\tYfin=%f\n",Yfin);
            MRh[i][ic]=Yfin/division;
        }
    }
    if(ic==10 || (ic ==11 && iPlot11==1)){
        // PLOT Rh /T
        QwtPlot *RhPlot=new QwtPlot();
        QwtPlotCurve *curve1=new QwtPlotCurve("Curve 1");
        if(ic==10){
            RhPlot -> setTitle("R-hemispherical spectrum");
            RhPlot -> setAxisTitle(0,"Rhemispherical");
        }
        else{
            RhPlot -> setTitle("Transmittance spectrum");
            RhPlot -> setAxisTitle(0,"T");
        }
        RhPlot -> setAxisTitle(2,"wavelength (nm)");
        RhPlot -> setAxisScale(2,MRh[0][0],MRh[Nw-1][0],0.);
        for(int ii=0;ii< Nw;ii++){
            if(ic==11)
                MRh[ii][ic]=MRh[ii][ic]/division;
            xPlot[ii]=MRh[ii][0];
            yPlot[ii]=MRh[ii][ic];
        }
        curve1->setSamples(xPlot, yPlot, Nw);
        curve1->attach(RhPlot);
        // Make the grid on
        QwtPlotGrid *grid = new QwtPlotGrid();
        grid->setPen(QPen(Qt::gray, 0.0, Qt::DotLine));
        grid->enableX(true);
        grid->enableXMin(true);
        grid->enableY(true);
        grid->enableYMin(true);
        grid->attach(RhPlot);
        RhPlot->replot();
        RhPlot->show();
    }
}


void nkSUB(double lambda, double T, double R){
    double Rf=(2.+T*T-pow(1.-R,2.)-
               sqrt( pow(2.+T*T-pow(1.-R,2.),2.)-4.*R*(2.-R)))/
               (2.*(2.-R));
    kOl=lambda/(4.*Pig*dW)*log(Rf*T/(R-Rf));
    nOl=(1.+Rf)/(1.-Rf);
    nOl=nOl+sqrt(4.*Rf/pow(1.-Rf,2.)-kOl*kOl);//"+" choice
    //printf("nkSUB: T=%f R=%f Rf=%f ArgSqrt=%f n=%f k=%e\n",T,R,Rf,4.*Rf/pow(1.-Rf,2.)-kOl*kOl,nOl,kOl);
    //fflush(stdout);
}



void RsRpEMA(int iWL, double theta){
    double alpha,T1s,T1p,R1s,R1p;
    complex<double> theta1;
    complex<double> N0(1.0,0.0);//air
    complex<double> theta0(theta,0.0);
    complex<double> N1(MRh[iWL][6],-MRh[iWL][7]);//nk-window
    RsRp(N0,N1,theta0,T1s,T1p,R1s,R1p,theta1);
    alpha=4.0*Pig/MRh[iWL][0]*imag(-sqrt(N1*N1-N0*N0*sin(theta0)*sin(theta0)));
    TEMAs=T1s*T1s/(exp(2.*alpha*dW)-R1s*R1s);
    TEMAp=T1p*T1p/(exp(2.*alpha*dW)-R1p*R1p);
    TEMAu=0.5*(TEMAs+TEMAp);
    REMAs=R1s+pow(1.0-R1s,2.0)*R1s/(exp(2.*alpha*dW)-R1s*R1s);
    REMAp=R1p+pow(1.0-R1p,2.0)*R1p/(exp(2.*alpha*dW)-R1p*R1p);
    REMAu=0.5*(REMAs+REMAp);
}


void RsRp(complex<double> N1, complex<double> N2, complex<double> theta1,
          double& Ts,double& Tp,double& Rs,double& Rp, complex<double>& theta2){
    theta2=asin(N1*sin(theta1)/N2);
    Ts=4.*real(N1*cos(theta1))*real(N2*cos(theta2))/pow(abs(N1*cos(theta1)+N2*cos(theta2)),2.);
    Tp=4.*real(N2*cos(theta1))*real(N1*cos(theta2))/pow(abs(N2*cos(theta1)+N1*cos(theta2)),2.);
    Rs=pow(abs((N1*cos(theta1)-N2*cos(theta2))/(N1*cos(theta1)+N2*cos(theta2))),2.);
    Rp=pow(abs((N2*cos(theta1)-N1*cos(theta2))/(N2*cos(theta1)+N1*cos(theta2))),2.);
    //printf("RsRp: N1=%f-i%f N2=%f-i%f Ts=%f Tp=%f Rs=%f Rp=%f\n",
    //       real(N1),imag(N1),real(N2),imag(N2),Ts,Tp,Rs,Rp);
    //fflush(stdout);
    //  if(Rs<0. || Rp<0. || Rs>1.0 || Rp>1.0){
    //   cout << "N1= " << N1 << "\n";
    //   cout << "N2= " << N2 << "\n";
    //   cout << "cos(theta1)= " << cos(theta1) << "\n";
    //   cout << "cos(theta2)= " << cos(theta2) << "\n";
    //  }
}
