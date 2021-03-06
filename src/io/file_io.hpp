#ifndef RAXML_FILE_IO_HPP_
#define RAXML_FILE_IO_HPP_

#include <fstream>

#include "../Tree.hpp"
#include "../AncestralStates.hpp"
#include "../bootstrap/BootstrapTree.hpp"
#include "../bootstrap/TransferBootstrapTree.hpp"
#include "../bootstrap/BootstrapGenerator.hpp"
#include "../PartitionedMSAView.hpp"

class NewickStream : public std::fstream
{
public:
  NewickStream(const std::string& fname) : std::fstream(fname, std::ios::out), _brlens(true) {};
  NewickStream(const std::string& fname, std::ios_base::openmode mode) :
    std::fstream(fname, mode), _brlens(true) {};

  bool brlens() const { return _brlens; }
  void brlens(bool v) { _brlens = v; }

private:
  bool _brlens;
};

class TBEExtraTableStream : public std::fstream
{
public:
	TBEExtraTableStream(const std::string& fname) : std::fstream(fname, std::ios::out) {};
	TBEExtraTableStream(const std::string& fname, std::ios_base::openmode mode) :
    std::fstream(fname, mode) {};
};

class TBEExtraArrayStream : public std::fstream
{
public:
	TBEExtraArrayStream(const std::string& fname) : std::fstream(fname, std::ios::out) {};
	TBEExtraArrayStream(const std::string& fname, std::ios_base::openmode mode) :
    std::fstream(fname, mode) {};
};

class TBEExtraTreeStream : public std::fstream
{
public:
	TBEExtraTreeStream(const std::string& fname) : std::fstream(fname, std::ios::out) {};
	TBEExtraTreeStream(const std::string& fname, std::ios_base::openmode mode) :
    std::fstream(fname, mode) {};
};

class MSAFileStream
{
public:
  MSAFileStream(const std::string& fname) :
    _fname(fname) {}

  const std::string& fname() const { return _fname; };

private:
  std::string _fname;
};

class PhylipStream : public MSAFileStream
{
public:
  PhylipStream(const std::string& fname, bool interleaved = true) :
    MSAFileStream(fname), _interleaved(interleaved) {}

  bool interleaved() const { return _interleaved; }
private:
  bool _interleaved;
};

class FastaStream : public MSAFileStream
{
public:
  FastaStream(const std::string& fname) : MSAFileStream(fname) {}
};

class CATGStream : public MSAFileStream
{
public:
  CATGStream(const std::string& fname) : MSAFileStream(fname) {}
};

class RBAStream : public MSAFileStream
{
public:
  RBAStream(const std::string& fname) : MSAFileStream(fname) {}

  static bool rba_file(const std::string& fname, bool check_version = false);
};

class RaxmlPartitionStream : public std::fstream
{
public:
  RaxmlPartitionStream(const std::string& fname, bool use_range_string = false) :
    std::fstream(fname, std::ios::out), _offset(0), _print_model_params(false),
    _use_range_string(use_range_string) {}
  RaxmlPartitionStream(const std::string& fname, std::ios_base::openmode mode) :
    std::fstream(fname, mode), _offset(0), _print_model_params(false), _use_range_string(false) {}

  bool print_model_params() const { return _print_model_params; }
  void print_model_params(bool value) { _print_model_params = value; }

  void reset() { _offset = 0; }
  void put_range(const PartitionInfo& part_info)
  {
    if (_use_range_string && !part_info.range_string().empty())
      *this << part_info.range_string();
    else
    {
      auto part_len = part_info.msa().num_sites();
      *this << (_offset+1) << "-" << (_offset+part_len);
      _offset += part_len;
    }
  }

private:
  size_t _offset;
  bool _print_model_params;
  bool _use_range_string;
};

class FileIOStream : public std::fstream
{
public:
  FileIOStream(const std::string& fname, std::ios_base::openmode mode = std::ios::out) :
    std::fstream(fname, mode), _delim("\t"), _precision(6) {};

  const std::string& delim() { return _delim; };
  void delim(const std::string& del) { _delim = del; };
  unsigned int precision() { return _precision; };
  void precision(unsigned int prec) { _precision = prec; };

protected:
  std::string _delim;
  unsigned int _precision;
};

class AncestralProbStream : public FileIOStream
{
public:
  AncestralProbStream(const std::string& fname) : FileIOStream(fname) {};
  AncestralProbStream(const std::string& fname, std::ios_base::openmode mode) :
    FileIOStream(fname, mode) {};
};

class AncestralStateStream : public FileIOStream
{
public:
  AncestralStateStream(const std::string& fname) : FileIOStream(fname) {};
  AncestralStateStream(const std::string& fname, std::ios_base::openmode mode) :
    FileIOStream(fname, mode) {};
};

NewickStream& operator<<(NewickStream& stream, const pll_unode_t& root);
NewickStream& operator<<(NewickStream& stream, const pll_utree_t& tree);
NewickStream& operator<<(NewickStream& stream, const Tree& tree);
NewickStream& operator>>(NewickStream& stream, Tree& tree);

NewickStream& operator<<(NewickStream& stream, const AncestralStates& ancestral);
//NewickStream& operator>>(NewickStream& stream, BootstrapTree& tree);

NewickStream& operator<<(NewickStream& stream, const SupportTree& tree);

PhylipStream& operator>>(PhylipStream& stream, MSA& msa);
FastaStream& operator>>(FastaStream& stream, MSA& msa);
CATGStream& operator>>(CATGStream& stream, MSA& msa);
MSA msa_load_from_file(const std::string &filename, const FileFormat format);

PhylipStream& operator<<(PhylipStream& stream, const MSA& msa);
PhylipStream& operator<<(PhylipStream& stream, const PartitionedMSA& msa);
PhylipStream& operator<<(PhylipStream& stream, const PartitionedMSAView& msa);
PhylipStream& operator<<(PhylipStream& stream, const BootstrapMSA& bs_msa);

RBAStream& operator<<(RBAStream& stream, const PartitionedMSA& part_msa);
RBAStream& operator>>(RBAStream& stream, PartitionedMSA& part_msa);

RaxmlPartitionStream& operator>>(RaxmlPartitionStream& stream, PartitionInfo& part_info);
RaxmlPartitionStream& operator>>(RaxmlPartitionStream& stream, PartitionedMSA& parted_msa);

RaxmlPartitionStream& operator<<(RaxmlPartitionStream& stream, const PartitionInfo& part_info);
RaxmlPartitionStream& operator<<(RaxmlPartitionStream& stream, const PartitionedMSA& parted_msa);
RaxmlPartitionStream& operator<<(RaxmlPartitionStream& stream, const PartitionedMSAView& parted_msa);

AncestralProbStream& operator<<(AncestralProbStream& stream, const AncestralStates& ancestral);
AncestralStateStream& operator<<(AncestralStateStream& stream, const AncestralStates& ancestral);

TBEExtraTableStream& operator<<(TBEExtraTableStream& stream, const TransferBootstrapTree& tree);
TBEExtraArrayStream& operator<<(TBEExtraArrayStream& stream, const TransferBootstrapTree& tree);
TBEExtraTreeStream& operator<<(TBEExtraTreeStream& stream, const TransferBootstrapTree& tree);

std::string to_newick_string_rooted(const Tree& tree, double root_brlen = 0.0);

#endif /* RAXML_FILE_IO_HPP_ */
