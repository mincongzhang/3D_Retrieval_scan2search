
// CMeshRetrievalMFC.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "MeshRetrievalMFC.h"
#include "MeshRetrievalMFCDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMeshRetrievalMFCApp

BEGIN_MESSAGE_MAP(CMeshRetrievalMFCApp, CWinAppEx)
	ON_COMMAND(ID_HELP, CWinAppEx::OnContextHelp)
END_MESSAGE_MAP()

// CMeshRetrievalMFCApp construction

CMeshRetrievalMFCApp::CMeshRetrievalMFCApp()
{

	// Place all significant initialization in InitInstance
}


// The one and only CMeshRetrievalMFCApp object

CMeshRetrievalMFCApp theApp;


// CMeshRetrievalMFCApp initialization

BOOL CMeshRetrievalMFCApp::InitInstance()
{
	// InitCommonControlsEx() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);


	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinAppEx::InitInstance();

	AfxEnableControlContainer();

	// Standard initialization
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));

	CMeshRetrievalMFCDlg dlg;
	m_pMainWnd = &dlg;
	INT_PTR nResponse = dlg.DoModal();
	if (nResponse == IDOK)
	{

		//  dismissed with OK
	}
	else if (nResponse == IDCANCEL)
	{

		//  dismissed with Cancel
	}

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.
	return FALSE;
}
