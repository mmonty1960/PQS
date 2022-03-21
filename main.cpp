/*
Particle Qualification software
version for SolarPACES diffusion PQS

Author: Marco Montecchi
        ENEA-Casaccia
        marco.montecchi@enea.it

This file is part of PQS.

    PQSexpo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation version 3 of the License

    PQSexpo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2012-2022 Marco Montecchi
*/
#include <QApplication>
#include "PQS.h"
 
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    PQS w;
    w.show();
    return app.exec();
}

