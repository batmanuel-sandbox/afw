// -*- lsst-c++ -*-
#ifndef AFW_TABLE_Simple_h_INCLUDED
#define AFW_TABLE_Simple_h_INCLUDED



#include "lsst/afw/table/RecordInterface.h"
#include "lsst/afw/table/TableInterface.h"

namespace lsst { namespace afw { namespace table {

class SimpleRecord;
class SimpleTable;

/// @brief A tag class for SimpleTable and SimpleRecord to be used with the interface classes.
struct Simple {
    typedef SimpleRecord Record;
    typedef SimpleTable Table;
};

/**
 *  @brief A bare-bones record class intended for testing and generic tabular data.
 */
class SimpleRecord : public RecordInterface<Simple> {
public:

    //@{
    /**
     *  @brief Add a new child record to the same table this record belongs to.
     *
     *  Will throw LogicErrorException if !isLinked().
     */
    SimpleRecord addChild() const { return _addChild(); }
    SimpleRecord addChild(RecordId id) const { return _addChild(id); }
    //@}

private:

    friend class detail::Access;

    SimpleRecord(RecordBase const & other) : RecordInterface<Simple>(other) {}
};

/**
 *  @brief A bare-bones record class intended for testing and generic tabular data.
 */
class SimpleTable : public TableInterface<Simple> {
public:

    /**
     *  @brief Construct a new table.
     *
     *  @param[in] layout            Layout that defines the fields, offsets, and record size for the table.
     *  @param[in] capacity          Number of records to pre-allocate space for in the first block.  This
     *                               overrides nRecordsPerBlock for the first block and the first block only.
     *  @param[in] nRecordsPerBlock  Number of records to allocate space for in each block.  This is almost
     *                               entirely a performance-only consideration, but it does affect whether
     *                               a table will be remain consolidated after adding records.
     *  @param[in] idFactory         Factory class to generate record IDs when they are not explicitly given.
     *                               If empty, defaults to a simple counter that starts at 1.
     */
    SimpleTable(
        Layout const & layout,
        int capacity = 0,
        int nRecordsPerBlock = 256,
        IdFactory::Ptr const & idFactory = IdFactory::Ptr()
    ) : TableInterface<Simple>(layout, capacity, nRecordsPerBlock, idFactory) {}

    
    /// @brief Create and add a new record with an ID generated by the table's IdFactory.
    Record addRecord() const { return _addRecord(); }

    /// @brief Create and add a new record with an explicit RecordId.
    Record addRecord(RecordId id) const { return _addRecord(id); }

};

}}} // namespace lsst::afw::table

#endif // !AFW_TABLE_Simple_h_INCLUDED
