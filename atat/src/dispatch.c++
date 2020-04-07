#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <fstream.h>
#include "misc.h"
#include "array.h"
#include "stringo.h"
#include "mpiinterf.h"

main(int argc, char **argv)
{
  MyMPIobj.init(argc,argv);
  MyMPIobj.barrier();

  if (MyMPIobj.is_root()) {
    ifstream file(argv[1]);
    Array<int> done(MyMPIobj.numproc-1);
    zero_array(&done);
    int numdone=0;
    while (numdone<MyMPIobj.numproc-1) {
      int who;
      //cout << "check" << endl;
      MyMPI_RecvStream(&who,MPI_ANY_SOURCE,0);
      //cout << "got " << who << endl;
      char buf[MAX_LINE_LEN];
      char c;
      file.get(buf,MAX_LINE_LEN-1);
      file.get(c);
      AutoString str(buf);
      //cout << "read:" << strlen(buf) << ":" << buf << endl;
      MyMPI_SendStream(&str,who,who);
      if (strlen(buf)==0) {
	done(who-1)=1;
	numdone++;
      }
    }
  }
  else {
    while (1) {
      MyMPI_SendStream(&(MyMPIobj.id),0,0);
      //      cout << "sent to root" << MyMPIobj.id << endl;
      AutoString str;
      MyMPI_RecvStream(&str,0,MyMPIobj.id);
      if (str.len()==0) break;
      cout << "dispatch: Process " << MyMPIobj.id << " of " << MyMPIobj.numproc-1 << " is running " << str << endl;
      system(str);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "All done." << endl;
}
