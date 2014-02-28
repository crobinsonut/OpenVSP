import vsp

# class Geometry(object):
# 	def __init__(self, geometry_id):
# 		self._geometry_id = geometry_id

#     @property
# 	def name(self):
# 		vsp.GetGeomName(self._geometry_id)

# 	@property.setter
# 	def name(self, value):
# 		vsp.SetGeomName(self._geometry_id, value)

class Parameter(object):
	def __init__(self, parameter_id):
		self._parameter_id = parameter_id

	def __get__(self, instance, owner):
		return vsp.GetParmVal(self._parameter_id)

	def __set__(self, instance, value):
		return vsp.SetParmVal(self._parameter_id, value, True)

class ParameterGroup(object):
	def __init__(self, geometry_id, group_name, parameter_names):
		for parameter_name in parameter_names:
			parameter_id = vsp.GetParm(geometry_id, parameter_name, group_name)

			setattr(self.__class__, parameter_name.lower(), Parameter(parameter_id))

class Design(ParameterGroup):
	def __init__(self, geometry_id):
		super(Design, self).__init__(
											geometry_id,
											"Design",
											('Length', 'XSecConnect',),
									)

class BoundingBox(ParameterGroup):
	def __init__(self, geometry_id):
		super(BoundingBox, self).__init__(
											geometry_id,
											"BBox",
											('X_Min', 'X_len',
										     'Y_Min', 'Y_len',
										     'Z_Min', 'Z_len',)
											)

if __name__ == "__main__":
	pod_id = vsp.AddGeom("POD")
	design = Design(pod_id)
	print design.length
	print design.xsecconnect

	design.length = 5.0
	design.xsecconnect = 20.0

	print design.length
	print design.xsecconnect

