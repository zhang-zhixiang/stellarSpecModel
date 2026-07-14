class SpecModelError(Exception):
    """stellarSpecModel 基础异常类"""
    pass

class GridPointNotFoundError(SpecModelError):
    """恒星参数不在现有网格上时抛出"""
    pass

class ParameterConflictError(SpecModelError):
    """参数设置相互矛盾时抛出"""
    pass

class WavelengthOutOfBoundsError(SpecModelError):
    """重采样波长超出了模板波长物理边界时抛出"""
    pass

class AliasAlreadyExistsError(SpecModelError):
    """别名已经存在且未开启覆写时抛出"""
    pass