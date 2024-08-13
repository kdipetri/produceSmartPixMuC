import ROOT 
from math import sqrt

# for ROOT plots
#ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.LoadMacro("atlasrootstyle/AtlasStyle.C")
#ROOT.gROOT.LoadMacro("atlasrootstyle/AtlasUtils.C")
#ROOT.gROOT.LoadMacro("atlasrootstyle/AtlasLabels.C")
#ROOT.SetAtlasStyle()

class PlotHelper:

	def __init__(self):
		self.histos_1d = {}
		self.histos_2d = {}
		self.plotdir = "plots";
		self.doPdf = False
		self.doPng = True
		self.c1 = ROOT.TCanvas("c1","",800,800);
		self.c2 = ROOT.TCanvas("c2","",800,600);
		self.prettyCanvas()

	def prettyCanvas( self ):
		# for 1D histos
		self.c1.SetTopMargin(0.05)
		self.c1.SetBottomMargin(0.15)
		self.c1.SetLeftMargin(0.15)
		self.c1.SetRightMargin(0.05)
		# for 2D histos
		self.c2.SetTopMargin(0.05)
		self.c2.SetBottomMargin(0.15)
		self.c2.SetLeftMargin(0.2)
		self.c2.SetRightMargin(0.2)

	def plot1D( self, name, title, x, nbinsx, xmin, xmax, weight=1.) :

		h = self.histos_1d.get(name)

		if h is not None: # Hist exists, fill
			h.Fill(x, weight)
		else : # Create hist
			h = ROOT.TH1F(name, title, nbinsx, xmin, xmax)
			h.Sumw2()
			h.Fill(x, weight)
			self.histos_1d[name] = h

		return

	def plot2D( self, name, title, x, y, nbinsx, xmin, xmax, nbinsy, ymin, ymax, weight=1.) :

		h = self.histos_2d.get(name)

		if h is not None: # Hist exists, fill
			h.Fill(x, y, weight)
		else : # Create hist
			h = ROOT.TH2F(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
			h.Sumw2()
			h.Fill(x, y, weight)
			self.histos_2d[name] = h

		return

	def draw1D( self, h, logy=0, drawopt="hist", norm=False): 
		self.c1.cd();
		self.c1.Clear();

		h.Draw(drawopt);
		h.SetLineColor(ROOT.kBlue+1)
		h.SetLineWidth(2)
		h.SetMinimum(0)
	
		if norm : h.Scale(1.0/h.Integral(0,-1));

		# deal with overflow 
		h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX() + 1) + h.GetBinContent(h.GetNbinsX()) );
		h.SetBinContent(h.GetNbinsX() + 1, 0);
		h.SetBinError(h.GetNbinsX(), sqrt( h.GetBinError( h.GetNbinsX()+1 )**2 + h.GetBinError( h.GetNbinsX() )**2  ) );
		h.SetBinError(h.GetNbinsX() + 1, 0);

		# deal with underflow
		h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1));
		h.SetBinContent(0, 0);
		h.SetBinError(1, sqrt( h.GetBinError(0)**2 + h.GetBinError(1)**2  ) );
		h.SetBinError(0, 0);
		
		self.c1.SetLogy(logy);
		h.Write();
		if self.doPng : self.c1.Print("{}/{}.png".format(self.plotdir, h.GetName()) );
		if self.doPdf : self.c1.Print("{}/{}.pdf".format(self.plotdir, h.GetName()) );
		
		return 


	def draw2D( self, h, logz=0, drawopt="colz", norm=False): 
		self.c2.cd();
		self.c2.Clear();

		h.Draw(drawopt);
		h.GetZaxis().SetTitleOffset(1.5)
	
		if norm : h.Scale(1.0/h.Integral(0,-1,0,-1));
		
		self.c2.SetLogy(logz);
		h.Write();
		if self.doPng : self.c2.Print("{}/{}.png".format(self.plotdir, h.GetName()) );
		if self.doPdf : self.c2.Print("{}/{}.pdf".format(self.plotdir, h.GetName()) );
		
		return 

	def drawAll( self, drawopt1D="hist", drawopt2D="colz") :
		for i in self.histos_1d: 
			self.draw1D(self.histos_1d[i],drawopt=drawopt1D)
		for i in self.histos_2d: 
			self.draw2D(self.histos_2d[i],drawopt=drawopt2D)   
