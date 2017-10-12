#ifndef BOSON_H
#define BOSON_H
#include "itensor/mps/siteset.h"

namespace itensor {

class BosonSite;

using Boson = BasicSiteSet<BosonSite>;

class BosonSite
    {
    IQIndex s;
    public:

    BosonSite() { }

    BosonSite(IQIndex I) : s(I) { }

    BosonSite(int n, Args const& args = Args::global()){
      s = IQIndex{
        nameint("Boson ",n),
        Index(nameint("Emp ",n),1,Site),QN("Nb=",0),
        Index(nameint("Occ1 ",n),1,Site),QN("Nb=",1),
        Index(nameint("Occ2 ",n),1,Site),QN("Nb=",2),
        Index(nameint("Occ3 ",n),1,Site),QN("Nb=",3),
        Index(nameint("Occ4 ",n),1,Site),QN("Nb=",4),
        Index(nameint("Occ5 ",n),1,Site),QN("Nb=",5),
        Index(nameint("Occ6 ",n),1,Site),QN("Nb=",6),
        Index(nameint("Occ7 ",n),1,Site),QN("Nb=",7),
      };
    }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "0" || state == "Emp") {return s(1);}
        else
        if(state == "1" || state == "Occ1") {return s(2);}
        else
        if(state == "2" || state == "Occ2") {return s(3);}
        else
        if(state == "3" || state == "Occ3") {return s(4);}
        else
        if(state == "4" || state == "Occ4") {return s(5);}
        else
        if(state == "5" || state == "Occ5") {return s(6);}
        else
        if(state == "6" || state == "Occ6") {return s(7);}
        else
        if(state == "7" || state == "Occ7") {return s(8);}
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        IQIndexVal Emp(s(1)),
                   EmpP(sP(1)),
                   Occ1(s(2)),
                   Occ1P(sP(2)),
                   Occ2(s(3)),
                   Occ2P(sP(3)),
                   Occ3(s(4)),
                   Occ3P(sP(4)),
                   Occ4(s(5)),
                   Occ4P(sP(5)),
                   Occ5(s(6)),
                   Occ5P(sP(6)),
                   Occ6(s(7)),
                   Occ6P(sP(7)),
                   Occ7(s(8)),
                   Occ7P(sP(8));

        IQTensor Op(dag(s),sP);

        if(opname == "N")
            {
            Op.set(Emp,EmpP,0);
            Op.set(Occ1,Occ1P,1);
            Op.set(Occ2,Occ2P,2);
            Op.set(Occ3,Occ3P,3);
            Op.set(Occ4,Occ4P,4);
            Op.set(Occ5,Occ5P,5);
            Op.set(Occ6,Occ6P,6);
            Op.set(Occ7,Occ7P,7);
            }
        else
        if(opname == "N-1")
            {
            Op.set(Emp,EmpP,0);
            Op.set(Occ1,Occ1P,0);
            Op.set(Occ2,Occ2P,1);
            Op.set(Occ3,Occ3P,2);
            Op.set(Occ4,Occ4P,3);
            Op.set(Occ5,Occ5P,4);
            Op.set(Occ6,Occ6P,5);
            Op.set(Occ7,Occ7P,6);
            }
        else
        if(opname == "A")
            {
            Op.set(Occ1,EmpP,1);
            Op.set(Occ2,Occ1P,std::sqrt(2));
            Op.set(Occ3,Occ2P,std::sqrt(3));
            Op.set(Occ4,Occ3P,std::sqrt(4));
            Op.set(Occ5,Occ4P,std::sqrt(5));
            Op.set(Occ6,Occ5P,std::sqrt(6));
            Op.set(Occ7,Occ6P,std::sqrt(7));
            }
        else
        if(opname == "Adag")
            {
            Op.set(Emp,Occ1P,1);
            Op.set(Occ1,Occ2P,std::sqrt(2));
            Op.set(Occ2,Occ3P,std::sqrt(3));
            Op.set(Occ3,Occ4P,std::sqrt(4));
            Op.set(Occ4,Occ5P,std::sqrt(5));
            Op.set(Occ5,Occ6P,std::sqrt(6));
            Op.set(Occ6,Occ7P,std::sqrt(7));
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };


} //namespace itensor

#endif
