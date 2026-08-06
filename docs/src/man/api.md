```@meta
CurrentModule = XAM
DocTestSetup = quote
    using XAM
end
```

# API Reference

## SAM API

### Public

```@autodocs
Modules = [XAM.SAM]
private = false
```

### Internal
```@autodocs
Modules = [XAM.SAM]
public = false
```

## BAM API

### Public

```@autodocs
Modules = [XAM.BAM]
private = false
```

### Internal
```@autodocs
Modules = [XAM.BAM]
public = false
```

## Flags API
The Flags API provides predicates to test properties of the XAM records encoded in the SAM/BAM flags. They can be used for both SAM and BAM records.

```@docs
XAM.ispaired(record::XAMRecord)
XAM.isproperpair(record::XAMRecord)
XAM.isunmapped(record::XAMRecord)
XAM.ismapped(record::XAMRecord)
XAM.isnextunmapped(record::XAMRecord)
XAM.isnextmapped(record::XAMRecord)
XAM.isreversecomplemented(record::XAMRecord)
XAM.isforwardstrand(record::XAMRecord)
XAM.ispositivestrand(record::XAMRecord)
XAM.isreversestrand(record::XAMRecord)
XAM.isnegativestrand(record::XAMRecord)
XAM.isnextreversecomplemented(record::XAMRecord)
XAM.isfirstsegment(record::XAMRecord)
XAM.isread1(record::XAMRecord)
XAM.islastsegment(record::XAMRecord)
XAM.isread2(record::XAMRecord)
XAM.issecondaryalignment(record::XAMRecord)
XAM.isqcfail(record::XAMRecord)
XAM.isduplicate(record::XAMRecord)
XAM.issupplementaryalignment(record::XAMRecord)
XAM.isprimaryalignment(record::XAMRecord)
```

