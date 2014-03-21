import vsp
class Geometry(object):

    def __init__(self, geometry_id):

        self._id = geometry_id
        self.xsec_surfs = []
            
        parameter_ids = vsp.GetGeomParmIDs(self._id)

        for parameter_id in parameter_ids:
            parameter_name = vsp.GetParmName(parameter_id)

            setattr(self.__class__, parameter_name, Parameter(parameter_id))

        for index in range(vsp.GetNumXSecSurfs(self._id)):
            self.xsec_surfs.append(XSecSurf(vsp.GetXSecSurf(self._id, index)))

    def __str__(self):
        string = "{geometry_type}<name={geometry_name}, id={geometry_id}>"
        string = string.format(
                            geometry_type=self.__class__.__name__,
                            geometry_name=self.name,
                            geometry_id=self.id,
                            )
        return string

    @property
    def name(self):
        return vsp.GetGeomName(self._id)

    @name.setter
    def name(self, geometry_name):
        vsp.SetGeomName(self._id, geometry_name)

    @property
    def id(self):
        return self._id

    @classmethod
    def add_geometry(cls, geometry_type, parent=None):
        if not parent:
            geometry_id = vsp.AddGeom(geometry_type)
        else:
            geometry_id = vsp.AddGeom(geometry_type, parent.id)

        geometry = cls(geometry_id)

        return geometry

GeometryType = type("GeometryType", (), {geometry_type:geometry_type for geometry_type in vsp.GetGeomTypes()})

class Parameter(object):

	def __init__(self, parameter_id):
            self._id = parameter_id
            parameter_types = {
                                0 : float,
                                1 : int,
                                2 : bool,
                                3 : str,
                              }
            self._type = vsp.GetParmType(self._id)
            self._type = parameter_types[self._type]

	def __get__(self, instance, owner):
            return vsp.GetParmVal(self._id)

	def __set__(self, instance, value):
            if not isinstance(value, self._type):
                raise RuntimeError("Parameter must be of type {}".format(self._type))
                    
            vsp.SetParmVal(self._id, value, True)

class XSecSurf(object):
    def __init__(self, xsec_surf_id):
        self._id = xsec_surf_id
        self.xsecs = XSecList(self._id)

    @property
    def id(self):
        return self._id


class XSec(object):
    def __init__(self, xsec_id):
        self._id = xsec_id

    def __getattr__(self, name):
        self._add_attribute(name)
        return getattr(self, name)

    def __setattr__(self, name, value):
        if name != '_id':
            self._add_attribute(name)

        super(XSec, self).__setattr__(name, value)

    def _add_attribute(self, name):
        parameter_id = vsp.GetXSecParm(self._id, name)
        setattr(self.__class__, name, Parameter(parameter_id))

    @property
    def width(self):
        return vsp.GetXSecWidth(self._id)

    @width.setter
    def width(self, width):
        vsp.SetXSecWidthHeight(self._id, width, self.height)

    @property
    def height(self):
        return vsp.GetXSecHeight(self._id)

    @height.setter
    def height(self, height):
        vsp.SetXSecWidthHeight(self._id, height, self.width)

    @property
    def type(self):
        return vsp.GetXSecType(self._id)

    def set_points(self, point_vector):
        vsp.SetXSecPnts(self._id, point_vector)

class XSecList(object):
    def __init__(self, xsec_surf_id):
        self._xsec_surf_id = xsec_surf_id

    def __getitem__(self, key):
        return XSec(vsp.GetXSec(self._xsec_surf_id, key))
