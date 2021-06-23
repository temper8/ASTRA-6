     In Astra version 6.1 and later a new algorithm is implemented.
 Since that the electron-ion equipartition term can be treated implicitly.
     Although all previous commands remain valid it is recommended to
 switch to the new option where a request for including the electron-ion 
 equipartition is indicated by a star in the equation command lines.
 In this case all relevant terms are added by the Astra compiler and
 should not be explicitly included in a model.
     Examples:

 Old style. Explicit scheme
	TE:;	PE = PEICL + ...;	PET = ...
	TI:;	PI =-PEICL + ...;	PIT = ...
 Old style. Quasi-implicit scheme
	TE:;	PE = PEI*TI + ...;	PET = -PEI + ...
	TI:;	PI = PEI*TE + ...;	PIT = -PEI + ...

 New style. True-implicit scheme
	TE*:;	PE = ...;	PET = ...;
	TI*:;	PI = ...;	PIT = ...;

