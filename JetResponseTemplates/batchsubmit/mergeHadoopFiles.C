
void mergeHadoopFiles(const TString& indir, const TString& outpath) {
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");
  
  //gSystem->Load("../libMiniFWLite.so");
  TChain *chain = new TChain("Events");
  chain->SetMaxTreeSize(5000000000LL); //default is 100000000000LL

  std::cout << "Merging files from dir: " << indir << std::endl
	    << "Ouputting to: " << outpath << std::endl;

  chain->Add(indir + "/*.root");
  chain->Merge(outpath, "fast");

  std::cout << "Total events merged: " << chain->GetEntries() << std::endl;
  bmark->Stop("benchmark");
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f s", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f s", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
}
