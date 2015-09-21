#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QThread>

#include <QDebug>
#include <QDropEvent>
#include <QUrl>
#include <QList>
#include <QMimeData>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_pBufferCout(new ExtendedBuffer()),
    m_pBufferCerr(new ExtendedBuffer())
{
    ui->setupUi(this);

    m_pBufferCout->setTextColor(QColor( "green" ));
    m_pBufferCerr->setTextColor(QColor( "red" ));

    std::cout.rdbuf(m_pBufferCout);
    std::cerr.rdbuf(m_pBufferCerr);

    qRegisterMetaType<QTextCursor>("QTextCursor");
    qRegisterMetaType<QString>("QString");


    QObject::connect(m_pBufferCout, SIGNAL(sendText(QString,QColor)), this , SLOT(receiveText(QString,QColor)), Qt::QueuedConnection );
    QObject::connect(m_pBufferCerr, SIGNAL(sendText(QString,QColor)), this , SLOT(receiveText(QString,QColor)), Qt::QueuedConnection );

    QThread::currentThread()->setPriority(QThread::LowestPriority);

    m_pThread = new GTXthread( );
    QObject::connect(m_pThread, SIGNAL(getThreadFinished(void*)), this , SLOT( getApplicationFinished(void*)), Qt::QueuedConnection );

    this->setAcceptDrops(true);

}

MainWindow::~MainWindow()
{
    delete ui;

    delete m_pThread;
    //delete m_pBufferCerr;
    //delete m_pBufferCout;
}

void MainWindow::on_oButtonStartEC_clicked()
{

    QString oGTFpart = QString("-gtx ");
    oGTFpart = oGTFpart.append( m_sInputFile );

    QString sCL = oGTFpart;

    if (ui->oStatsBox->isChecked() == true)
    {
        QString oStatsPart = QString("-stats ");
        oStatsPart = oStatsPart.append(m_sStatsFile);

        sCL = sCL.append(" ").append(oStatsPart);
    }

    if (ui->oValidateBox->isChecked() == true)
    {
        sCL = sCL.append(" --validate");
    }


    m_pThread->setCL(sCL);

    m_pThread->start();
}

void MainWindow::receiveText(QString sString, QColor oColor)
{
    this->ui->textEdit->moveCursor(QTextCursor::End);

    this->ui->textEdit->setTextColor(oColor);
    this->ui->textEdit->insertPlainText(sString);

    qDebug() << sString;
}

void MainWindow::getApplicationFinished(void* pData)
{

    if (pData != NULL)
        free( pData );
}

void MainWindow::dropEvent(QDropEvent *oDropEvent)
{
    QList<QUrl> oUrls = oDropEvent->mimeData()->urls();
    foreach(QUrl oUrl, oUrls)
    {
        m_sInputFile = oUrl.path();

    }

    std::cout << "File Dropped: " << m_sInputFile.toStdString();

    ui->label_13->setText( m_sInputFile );

    return;

}

void MainWindow::dragEnterEvent(QDragEnterEvent *oDragEnterEvent)
{
    oDragEnterEvent->accept();
}

void MainWindow::on_oButtonConfig_clicked()
{
    m_sInputFile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                     "",
                                                     tr("Files (*.*)"));

        ui->label_13->setText( m_sInputFile );
}

void MainWindow::on_oButtonStats_clicked()
{
    m_sStatsFile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                     "",
                                                     tr("Files (*.*)"));

        ui->label_15->setText( m_sStatsFile );
}

void MainWindow::on_pushButton_clicked()
{
    if (!m_pThread->isFinished())
    {
        m_pThread->terminate();
        m_pThread->exit(0);

    }
}
