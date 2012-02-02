// -*- lsst-c++ -*-
#include <typeinfo>

#include "lsst/afw/table/io/FitsWriter.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/detail/Access.h"

#define SAVE_MEAS_SLOT(NAME, Name, TYPE, Type)                          \
    if (table->get ## Name ## Type ## Key().isValid()) {                \
        std::string s = table->getSchema().find(table->get ## Name ## Type ## Key()).field.getName(); \
        std::replace(s.begin(), s.end(), '.', '_');                     \
        _fits->writeKey(#NAME "_" #TYPE "_SLOT", s.c_str(), "Defines the " #Name #Type " slot"); \
    }                                                                   \
    if (table->get ## Name ## Type ## ErrKey().isValid()) {             \
        std::string s = table->getSchema().find(table->get ## Name ## Type ## ErrKey()).field.getName(); \
        std::replace(s.begin(), s.end(), '.', '_');                     \
        _fits->writeKey(#NAME "_" #TYPE "_ERR_SLOT", s.c_str(),         \
                        "Defines the " #Name #Type "Err slot");         \
    }                                                                   \
    if (table->get ## Name ## Type ## Flag ## Key().isValid()) {        \
        std::string s = table->getSchema().find(table->get ## Name ## Type ## FlagKey()).field.getName(); \
        std::replace(s.begin(), s.end(), '.', '_');                     \
        _fits->writeKey(#NAME "_" #TYPE "_FLAG_SLOT", s.c_str(),        \
                        "Defines the " #Name #Type "Flag slot");        \
    }

#define SAVE_FLUX_SLOT(NAME, Name) SAVE_MEAS_SLOT(NAME, Name, FLUX, Flux)
#define SAVE_CENTROID_SLOT() SAVE_MEAS_SLOT(, , CENTROID, Centroid)
#define SAVE_SHAPE_SLOT() SAVE_MEAS_SLOT(, , SHAPE, Shape)

#define LOAD_MEAS_SLOT(NAME, Name, TYPE, Type)                          \
    {                                                                   \
        std::string s, sErr, sFlag;                                     \
        _fits->readKey(#NAME "_" #TYPE "_SLOT", s);                     \
        _fits->readKey(#NAME "_" #TYPE "_ERR_SLOT", sErr);              \
        _fits->readKey(#NAME "_" #TYPE "_FLAG_SLOT", sFlag);            \
        if (_fits->status == 0) {                                       \
            std::replace(s.begin(), s.end(), '_', '.');                 \
            std::replace(sErr.begin(), sErr.end(), '_', '.');           \
            std::replace(sFlag.begin(), sFlag.end(), '_', '.');           \
            table->define ## Name ## Type(schema[s], schema[sErr], schema[sFlag]); \
        } else {                                                        \
            _fits->status = 0;                                          \
        }                                                               \
    }
    

#define LOAD_FLUX_SLOT(NAME, Name) LOAD_MEAS_SLOT(NAME, Name, FLUX, Flux)
#define LOAD_CENTROID_SLOT() LOAD_MEAS_SLOT(, , CENTROID, Centroid)
#define LOAD_SHAPE_SLOT() LOAD_MEAS_SLOT(, , SHAPE, Shape)

namespace lsst { namespace afw { namespace table {

namespace {

class SourceTableImpl;

class SourceRecordImpl : public SourceRecord {
public:

    explicit SourceRecordImpl(PTR(SourceTable) const & table) : SourceRecord(table) {}

};

class SourceTableImpl : public SourceTable {
public:

    explicit SourceTableImpl(
        Schema const & schema,
        PTR(IdFactory) const & idFactory
    ) : SourceTable(schema, idFactory) {}

    SourceTableImpl(SourceTableImpl const & other) : SourceTable(other) {}

private:

    virtual PTR(BaseTable) _clone() const {
        return boost::make_shared<SourceTableImpl>(*this);
    }

    virtual PTR(BaseRecord) _makeRecord() {
        PTR(SourceRecord) record = boost::make_shared<SourceRecordImpl>(getSelf<SourceTableImpl>());
        record->setId((*getIdFactory())());
        return record;
    }

};

class SourceFitsWriter : public io::FitsWriter {
public:

    explicit SourceFitsWriter(Fits * fits) : io::FitsWriter(fits) {}

protected:
    
    virtual void _writeTable(CONST_PTR(BaseTable) const & table);

    virtual void _writeRecord(BaseRecord const & record);

private:
    int _spanCol;
    int _peakCol;
};

class SourceFitsReader : public io::FitsReader {
public:

    explicit SourceFitsReader(Fits * fits) : io::FitsReader(fits), _spanCol(-1), _peakCol(-1) {}

protected:

    virtual Schema _readSchema(int nCols=-1);

    virtual PTR(BaseTable) _readTable(Schema const & schema);

    virtual PTR(BaseRecord) _readRecord(PTR(BaseTable) const & table);

private:
    int _spanCol;
    int _peakCol;
};

void SourceFitsWriter::_writeTable(CONST_PTR(BaseTable) const & t) {
    CONST_PTR(SourceTable) table = boost::dynamic_pointer_cast<SourceTable const>(t);
    if (!table) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot use a SourceFitsWriter on a non-Source table."
        );
    }
    io::FitsWriter::_writeTable(table);
    _spanCol = _fits->addColumn<int>("spans", 0, "footprint spans (y, x0, x1)");
    _peakCol = _fits->addColumn<float>("peaks", 0, "footprint peaks (fx, fy, peakValue)");
    _fits->writeKey("SPANCOL", _spanCol + 1, "Column with footprint spans.");
    _fits->writeKey("PEAKCOL", _peakCol + 1, "Column with footprint peaks (float values).");
    _fits->writeKey("AFW_TYPE", "SOURCE", "Tells lsst::afw to load this as a Source table.");
    SAVE_FLUX_SLOT(PSF, Psf);
    SAVE_FLUX_SLOT(MODEL, Model);
    SAVE_FLUX_SLOT(AP, Ap);
    SAVE_FLUX_SLOT(INST, Inst);
    SAVE_CENTROID_SLOT();
    SAVE_SHAPE_SLOT();
}

void SourceFitsWriter::_writeRecord(BaseRecord const & r) {
    SourceRecord const & record = static_cast<SourceRecord const &>(r);
    io::FitsWriter::_writeRecord(record);
    if (record.getFootprint()) {
        Footprint::SpanList const & spans = record.getFootprint()->getSpans();
        Footprint::PeakList const & peaks = record.getFootprint()->getPeaks();
        if (!spans.empty()) {
            std::vector<int> vec;
            vec.reserve(3 * spans.size());
            for (Footprint::SpanList::const_iterator j = spans.begin(); j != spans.end(); ++j) {
                vec.push_back((**j).getY());
                vec.push_back((**j).getX0());
                vec.push_back((**j).getX1());
            }
            _fits->writeTableArray(_row, _spanCol, vec.size(), &vec.front());
        }
        if (!peaks.empty()) {
            std::vector<float> vec;
            vec.reserve(3 * peaks.size());
            for (Footprint::PeakList::const_iterator j = peaks.begin(); j != peaks.end(); ++j) {
                vec.push_back((**j).getFx());
                vec.push_back((**j).getFy());
                vec.push_back((**j).getPeakValue());}
            _fits->writeTableArray(_row, _peakCol, vec.size(), &vec.front());
        }
    }
}

Schema SourceFitsReader::_readSchema(int nCols) {
    _fits->readKey("SPANCOL", _spanCol);
    if (_fits->status == 0) {
        --_spanCol;
    } else {
        _fits->status = 0;
        _spanCol = -1;
    }
    _fits->readKey("PEAKCOL", _peakCol);
    if (_fits->status == 0) {
        --_peakCol;
    } else {
        _fits->status = 0;
        _peakCol = -1;
    }
    int maxCol = std::min(_spanCol, _peakCol);
    return io::FitsReader::_readSchema(maxCol);
}

PTR(BaseTable) SourceFitsReader::_readTable(Schema const & schema) {
    PTR(SourceTable) table =  SourceTable::make(schema);
    LOAD_FLUX_SLOT(PSF, Psf);
    LOAD_FLUX_SLOT(MODEL, Model);
    LOAD_FLUX_SLOT(AP, Ap);
    LOAD_FLUX_SLOT(INST, Inst);
    LOAD_CENTROID_SLOT();
    LOAD_SHAPE_SLOT();
    return table;
}

PTR(BaseRecord) SourceFitsReader::_readRecord(PTR(BaseTable) const & table) {
    PTR(SourceRecord) record = boost::static_pointer_cast<SourceRecord>(io::FitsReader::_readRecord(table));
    if (!record) return record;
    boost::static_pointer_cast<SourceTable>(table)->getIdFactory()->notify(record->getId());
    int spanElementCount = (_spanCol >= 0) ? _fits->getTableArraySize(_row, _spanCol) : 0;
    int peakElementCount = (_peakCol >= 0) ? _fits->getTableArraySize(_row, _peakCol) : 0;
    if (spanElementCount || peakElementCount) {
        PTR(Footprint) fp = boost::make_shared<Footprint>();
        if (spanElementCount) {
            if (spanElementCount % 3) {
                throw LSST_EXCEPT(
                    afw::fits::FitsError,
                    afw::fits::makeErrorMessage(
                        _fits->fptr, _fits->status,
                        boost::format("Number of span elements (%d) must divisible by 3 (row %d)")
                        % spanElementCount % _row
                    )
                );
            }
            std::vector<int> spanElements(spanElementCount);
            _fits->readTableArray(_row, _spanCol, spanElementCount, &spanElements.front());
            std::vector<int>::iterator j = spanElements.begin();
            while (j != spanElements.end()) {
                int y = *j++;
                int x0 = *j++;
                int x1 = *j++;
                fp->addSpan(y, x0, x1);
            }
        }
        if (peakElementCount) {
            if (peakElementCount % 3) {
                throw LSST_EXCEPT(
                    afw::fits::FitsError,
                    afw::fits::makeErrorMessage(
                        _fits->fptr, _fits->status,
                        boost::format("Number of peak elements (%d) must divisible by 3 (row %d)")
                        % peakElementCount % _row
                    )
                );
            }
            std::vector<float> peakElements(peakElementCount);
            _fits->readTableArray(_row, _peakCol, peakElementCount, &peakElements.front());
            std::vector<float>::iterator j = peakElements.begin();
            while (j != peakElements.end()) {
                float x = *j++;
                float y = *j++;
                float value = *j++;
                fp->getPeaks().push_back(boost::make_shared<detection::Peak>(x, y, value));
            }
        }
        record->setFootprint(fp);
    }
    return record;
}

static io::FitsReader::FactoryT<SourceFitsReader> sourceFitsReaderFactory("SOURCE");

} // anonymous

SourceRecord::SourceRecord(PTR(SourceTable) const & table) : BaseRecord(table) {}

void SourceRecord::_assign(BaseRecord const & other) {
    try {
        SourceRecord const & s = dynamic_cast<SourceRecord const &>(other);
        _footprint = s._footprint;
    } catch (std::bad_cast&) {}
}

PTR(SourceTable) SourceTable::make(
    Schema const & schema,
    PTR(IdFactory) const & idFactory
) {
    if (!checkSchema(schema)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Schema for Source must contain at least the keys defined by getMinimalSchema()."
        );
    }
    return boost::make_shared<SourceTableImpl>(schema, idFactory);
}

SourceTable::SourceTable(
    Schema const & schema,
    PTR(IdFactory) const & idFactory
) : BaseTable(schema), _idFactory(idFactory)
{
    if (!_idFactory) _idFactory = IdFactory::makeSimple();
}

SourceTable::SourceTable(SourceTable const & other) :
    BaseTable(other),
    _idFactory(other._idFactory->clone()),
    _slotFlux(other._slotFlux), _slotCentroid(other._slotCentroid), _slotShape(other._slotShape)
{}

SourceTable::MinimalSchema::MinimalSchema() {
    detail::Access::markPersistent(schema);
    id = schema.addField<RecordId>("id", "unique ID for source");
    parent = schema.addField<RecordId>("parent", "unique ID of parent source");
    sky = schema.addField<float>("sky", "sky background at location of source", "DN/pix");
    skyErr = schema.addField<float>("sky.err", "sky background uncertainty at location of source",
                                    "DN/pix");
    coord = schema.addField<Coord>("coord", "position of source in ra/dec", "radians");
}

SourceTable::MinimalSchema & SourceTable::getMinimalSchema() {
    static MinimalSchema it;
    return it;
}

PTR(io::FitsWriter) SourceTable::makeFitsWriter(io::FitsWriter::Fits * fits) const {
    return boost::make_shared<SourceFitsWriter>(fits);
}

KeyTuple<Centroid> addCentroidFields(
    Schema & schema,
    std::string const & name,
    std::string const & doc
) {
    KeyTuple<Centroid> keys;
    keys.meas = schema.addField<Centroid::MeasTag>(name, doc, "pixels");
    keys.err = schema.addField<Centroid::ErrTag>(
        name + ".err", "covariance matrix for " + name, "pixels^2"
    );
    keys.flag = schema.addField<Flag>(
        name + ".flags", "set if the " + name + " measurement did not fully succeed"
    );
    return keys;
}

KeyTuple<Shape> addShapeFields(
    Schema & schema,
    std::string const & name,
    std::string const & doc
) {
    KeyTuple<Shape> keys;
    keys.meas = schema.addField<Shape::MeasTag>(
        name, doc, "pixels^2"
    );
    keys.err = schema.addField<Shape::ErrTag>(
        name + ".err", "covariance matrix for " + name, "pixels^4"
    );
    keys.flag = schema.addField<Flag>(
        name + ".flags", "set if the " + name + " measurement failed"
    );
    return keys;
}

KeyTuple<Flux> addFluxFields(
    Schema & schema,
    std::string const & name,
    std::string const & doc
) {
    KeyTuple<Flux> keys;
    keys.meas = schema.addField<Flux::MeasTag>(
        name, doc, "dn"
    );
    keys.err = schema.addField<Flux::ErrTag>(
        name + ".err", "uncertainty for " + name, "dn"
    );
    keys.flag = schema.addField<Flag>(
        name + ".flags", "set if the " + name + " measurement failed"
    );
    return keys;
}

template class VectorT<SourceRecord>;
template class VectorT<SourceRecord const>;

}}} // namespace lsst::afw::table
