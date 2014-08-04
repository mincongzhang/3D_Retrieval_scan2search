// MeshRetrievalMFCDlg.cpp : implementation file
//

#include "stdafx.h"
#include "MeshRetrievalMFC.h"
#include "MeshRetrievalMFCDlg.h"



#ifdef _DEBUG
#define new DEBUG_NEW
#endif

int base = 0;

// CAboutDlg dialog used for App About
char fileName[256]={0};

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

	// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CMeshRetrievalMFCDlg dialog
CMeshRetrievalMFCDlg::CMeshRetrievalMFCDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMeshRetrievalMFCDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMeshRetrievalMFCDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CMeshRetrievalMFCDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_WM_SIZE()
	ON_WM_TIMER()
	ON_BN_CLICKED(Load, &CMeshRetrievalMFCDlg::OnBnClickedLoad)
	ON_BN_CLICKED(Noise, &CMeshRetrievalMFCDlg::OnBnClickedNoise)
	ON_BN_CLICKED(Normalize, &CMeshRetrievalMFCDlg::OnBnClickedNormalize)
	ON_BN_CLICKED(Spharm, &CMeshRetrievalMFCDlg::OnBnClickedSpharm)
	ON_BN_CLICKED(CANDIDATE1, &CMeshRetrievalMFCDlg::OnBnClickedCandidate1)
	ON_BN_CLICKED(CANDIDATE2, &CMeshRetrievalMFCDlg::OnBnClickedCandidate2)
	ON_BN_CLICKED(CANDIDATE3, &CMeshRetrievalMFCDlg::OnBnClickedCandidate3)
	ON_BN_CLICKED(CANDIDATE4, &CMeshRetrievalMFCDlg::OnBnClickedCandidate4)
	ON_BN_CLICKED(CANDIDATE5, &CMeshRetrievalMFCDlg::OnBnClickedCandidate5)
	ON_BN_CLICKED(CANDIDATE6, &CMeshRetrievalMFCDlg::OnBnClickedCandidate6)

	ON_BN_CLICKED(Rasterize, &CMeshRetrievalMFCDlg::OnBnClickedRasterize)
	ON_BN_CLICKED(BatchTransform, &CMeshRetrievalMFCDlg::OnBnClickedBatchtransform)
	ON_BN_CLICKED(Denoise, &CMeshRetrievalMFCDlg::OnBnClickedDenoise)
END_MESSAGE_MAP()


// CMeshRetrievalMFCDlg message handlers
BOOL CMeshRetrievalMFCDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	CRect rect;

	// Get size and position of the template textfield we created before in the dialog editor
	GetDlgItem(IDC_PIC_OPENGL)->GetWindowRect(rect);

	// Convert screen coordinates to client coordinates
	ScreenToClient(rect);

	// Create OpenGL Control window
	m_oglWindow.oglCreate(rect, this);

	// Setup the OpenGL Window's timer to render
	m_oglWindow.m_unpTimer = m_oglWindow.SetTimer(1, 1, 0);



	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CMeshRetrievalMFCDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CMeshRetrievalMFCDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CMeshRetrievalMFCDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


void CMeshRetrievalMFCDlg::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);

	// TODO: Add your message handler code here
	switch (nType)
	{
	case SIZE_RESTORED:
		{
			if (m_oglWindow.m_bIsMaximized)
			{
				m_oglWindow.OnSize(nType, cx, cy);
				m_oglWindow.m_bIsMaximized = false;
			}

			break;
		}

	case SIZE_MAXIMIZED:
		{
			m_oglWindow.OnSize(nType, cx, cy);
			m_oglWindow.m_bIsMaximized = true;

			break;
		}
	}
}

void CMeshRetrievalMFCDlg::OnTimer(UINT_PTR nIDEvent)
{
	// TODO: Add your message handler code here and/or call default
	CDialog::OnTimer(nIDEvent);
}

void CMeshRetrievalMFCDlg::OnBnClickedLoad()
{
	CFileDialog dlg(true,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		"all files(*.*)|*.*|model files(*.ply)|*.ply", NULL);
	if (dlg.DoModal()!= IDOK) return;

	MyMesh mesh;

	/*char fileName[256]={0};*/
	sprintf_s(fileName,"%s\0",dlg.GetPathName());

	try
	{
		if ( !OpenMesh::IO::read_mesh(mesh, fileName) )
		{
			std::cerr << "Cannot read mesh from the directory" << std::endl;
		}
	}
	catch( std::exception& x )
	{
		std::cerr << x.what() << std::endl;
	}
	meshQueue.reserve(10);
	meshQueue.push_back(mesh);
}

void CMeshRetrievalMFCDlg::OnBnClickedNoise()
{
	NOISE_CONTROL = true;
}

void CMeshRetrievalMFCDlg::OnBnClickedNormalize()
{
	NORMALIZE_CONTROL = true;
}

void CMeshRetrievalMFCDlg::OnBnClickedRasterize()
{
	RASTERIZE_CONTROL = true;
}

void CMeshRetrievalMFCDlg::OnBnClickedSpharm()
{
	SPHARM_CONTROL = true;
}

void CMeshRetrievalMFCDlg::OnBnClickedCandidate1()
{
	ChooseCandidate(candidate_index_array,base+0);
}

void CMeshRetrievalMFCDlg::OnBnClickedCandidate2()
{
	ChooseCandidate(candidate_index_array,base+1);
}

void CMeshRetrievalMFCDlg::OnBnClickedCandidate3()
{
	ChooseCandidate(candidate_index_array,base+2);
}

void CMeshRetrievalMFCDlg::OnBnClickedCandidate4()
{
	ChooseCandidate(candidate_index_array,base+3);
}

void CMeshRetrievalMFCDlg::OnBnClickedCandidate5()
{
	ChooseCandidate(candidate_index_array,base+4);
}

void CMeshRetrievalMFCDlg::OnBnClickedCandidate6()
{
	ChooseCandidate(candidate_index_array,base+5);
}



void CMeshRetrievalMFCDlg::OnBnClickedBatchtransform()
{
	BATCH_CONTROL = true;
}


void CMeshRetrievalMFCDlg::OnBnClickedDenoise()
{
	LAPLACE_DENOISE_CONTROL = true;
}
