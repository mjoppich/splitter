#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QThread>
#include <libPTA/PTA_ErrorCorrect.h>
#include <QDebug>
#include <QDropEvent>
#include <QUrl>
#include <QList>
#include <QMimeData>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_pBufferCout(new PTAExtendedBuffer()),
    m_pBufferCerr(new PTAExtendedBuffer())
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

    m_pThread = new PTAErrorCorrectionThread();
    QObject::connect(m_pThread, SIGNAL(getPTAECFinished(PTA::PTAStatistics*)), this , SLOT(getPTAECFinished(PTA::PTAStatistics*)), Qt::QueuedConnection );

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

    // first check whether there's a config file loaded!
    if (m_sConfigFile.compare("") == 0)
    {
        // empty config !

        std::cout << "ERROR: No Config File Set" << std::endl;
        m_sConfigFile = QString("/home/mjoppich/master/trunk/src/configs/dus12_16.config");

       // return;
    }

    int iConfigs = 1;
    const char** pConfigs = (const char**) malloc(sizeof(char*) * iConfigs);
    pConfigs[0] = m_sConfigFile.toStdString().c_str();

    m_pThread->setConfigs(pConfigs, iConfigs);
    m_pThread->start();
}

void MainWindow::receiveText(QString sString, QColor oColor)
{
    this->ui->textEdit->moveCursor(QTextCursor::End);

    this->ui->textEdit->setTextColor(oColor);
    this->ui->textEdit->insertPlainText(sString);

    qDebug() << sString;
}

void MainWindow::getPTAECFinished(PTA::PTAStatistics* pStats)
{



    this->ui->label_2->setText( QString::number(pStats->dSpentTime) );
    this->ui->label_8->setText( QString::number(pStats->dSpentTimeCorrection) );
    this->ui->label_9->setText( QString::number(pStats->iCorrectedReads) );
    this->ui->label_10->setText( QString::number(pStats->iTrimmedReads) );
    this->ui->label_11->setText( QString::number(pStats->iReadsInGraph) );
    this->ui->label_12->setText( QString::number(pStats->iSmersInGraph) );
    this->ui->label_22->setText( QString::number(pStats->iRoutesInGraph) );
    this->ui->label_24->setText( QString::number(pStats->iWidthInGraph) );


    delete pStats;
}

void MainWindow::dropEvent(QDropEvent *oDropEvent)
{
    QList<QUrl> oUrls = oDropEvent->mimeData()->urls();
    foreach(QUrl oUrl, oUrls)
    {
        m_sConfigFile = oUrl.path();

    }

    std::cout << "File Dropped: " << m_sConfigFile.toStdString();

    ui->label_13->setText( m_sConfigFile );

    return;

}

void MainWindow::dragEnterEvent(QDragEnterEvent *oDragEnterEvent)
{
    oDragEnterEvent->accept();
}

void MainWindow::on_oButtonConfig_clicked()
{
    m_sConfigFile = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                     "",
                                                     tr("Files (*.*)"));

    // if no file is loaded, we choose one by default for demonstration purposes
        if (m_sConfigFile.compare("") == 0)
    {
        // empty config !
        std::cerr << "ERROR: No Config File Set" << std::endl;
        std::cout << "ERROR: Config Set to dus12_16" << std::endl;

        m_sConfigFile = QString("/home/mjoppich/master/trunk/src/configs/dus12_16.config");
    }

        ui->label_13->setText( m_sConfigFile );
}

void MainWindow::on_pushButton_clicked()
{
    if (!m_pThread->isFinished())
    {
        m_pThread->terminate();
        m_pThread->exit(0);

    }
}
