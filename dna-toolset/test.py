#                                                      Classes

class Employee:
    def __init__(self, first, last, pay):
        self.first = first
        self.last = last
        self.pay = pay
        self.email = first + '.' + last + '@gmail.com'

    def fullname(self):
        # make sure to use (self): as the instance, or else you will get a type error
        return '{} {}'.format(self.first, self.last)
        # this is an attribute, therefore you need () when using print(emp_1.fullname())

emp_1 = Employee('Nick', 'Mabini', 60000)
emp_2 = Employee('Gerald', 'Jones', 50000)
emp_3 = Employee('Tobi', 'John', 40000)
emp_4 = Employee('Bill', 'Gate', 30000)

print(emp_1.email)
print(emp_2.email)
print(emp_3.email)
print(emp_4.email)

print(emp_1.fullname())
print(emp_2.fullname())
