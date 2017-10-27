#ifndef BOSON_H
#define BOSON_H
#include "itensor/mps/siteset.h"

namespace itensor {

class BosonSite;

template<typename SiteType>
class BosonSiteSet : public SiteSet
    {
    public:

    BosonSiteSet() { }

    BosonSiteSet(int N, int d)
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            sites.set(j,SiteType(j,d));
            }
        SiteSet::init(std::move(sites));
        }

    BosonSiteSet(std::vector<IQIndex> const& inds)
        {
        int N = inds.size();
        auto sites = SiteStore(N);
        for(int j = 1, i = 0; j <= N; ++j, ++i)
            {
            auto& Ii = inds.at(i);
            sites.set(j,SiteType(Ii));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        SiteSet::readType<SiteType>(s);
        }

  };

using Boson = BosonSiteSet<BosonSite>;

class BosonSite {
  IQIndex s;
  int d;
  std::vector<std::string> stateNames;

  public:

  BosonSite() { }

  BosonSite(IQIndex I) : s(I) { }

  BosonSite(int n, int d)
   : d(d) {
      char space[1] = { ' ' };
      auto v = stdx::reserve_vector<IndexQN>(1+d);
      v.emplace_back(Index(nameint("Emp ",n),1,Site),QN("Nb=",0));
      stateNames.emplace_back("Emp");

      for (size_t i = 1; i <= d; i++) {
        auto name = nameint("Occ",i);
        stateNames.emplace_back(name);
        name += space;
        v.emplace_back(Index(nameint(name,n),1,Site),QN("Nb=",i));
      }

      s = IQIndex(nameint("Boson ",n), std::move(v) );
    }

  IQIndex index() const { return s; }

  IQIndexVal state(std::string const& state) {
      for (size_t j = 0; j < stateNames.size(); j++) {
        if (state == stateNames.at(j) || state == nameint("",j)) {return s(j+1);}
      }

      Error("State " + state + " not recognized");
      return IQIndexVal{};
    }

	IQTensor op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        std::vector<IQIndexVal> indices(d+1);
        std::vector<IQIndexVal> indicesP(d+1);
        for (size_t j = 0; j <= d; j++) {
          indices.at(j) = s(j+1);
          indicesP.at(j) = sP(j+1);
        }

        IQTensor Op(dag(s),sP);

        if(opname == "N")
            {
              for (size_t j = 0; j <= d; j++) {
                  Op.set(indices.at(j),indicesP.at(j),j);
              }
            }
        else
        if(opname == "N-1")
            {
              for (size_t j = 1; j <= d; j++) {
                  Op.set(indices.at(j),indicesP.at(j),j-1);
              }
            }
        else
        if(opname == "A")
            {
              for (size_t j = 1; j <= d; j++) {
                  Op.set(indices.at(j),indicesP.at(j-1),std::sqrt(j));
              }
            }
        else
        if(opname == "Adag")
            {
              for (size_t j = 1; j <= d; j++) {
                  Op.set(indices.at(j-1),indicesP.at(j),std::sqrt(j));
              }
            }
        else
	if(opname == "N(N-1)")
	    {
	      for (size_t j = 1; j<= d; j++) {
		  Op.set(indices.at(j),indicesP.at(j),j*j-j);
	      }
	    }
        else
	if(opname == "I")
	    {
	      for (size_t j = 1; j<= d; j++) {
		  Op.set(indices.at(j),indicesP.at(j),1);
	      }
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
