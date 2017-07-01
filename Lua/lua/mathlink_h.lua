--[[
local f = io.open('mathlink.lua.h')
local r = f:read('*all')
f:close()
require'ffi'.cdef(r)
]]
--mathlink.h (defines made enums, macros removed)
require'ffi'.cdef[[

enum PKT {
ILLEGALPKT = 0,
CALLPKT = 7,
EVALUATEPKT = 13,
RETURNPKT = 3,
INPUTNAMEPKT = 8,
ENTERTEXTPKT = 14,
ENTEREXPRPKT = 15,
OUTPUTNAMEPKT = 9,
RETURNTEXTPKT = 4,
RETURNEXPRPKT = 16,
DISPLAYPKT = 11,
DISPLAYENDPKT = 12,
MESSAGEPKT = 5,
TEXTPKT = 2,
INPUTPKT = 1,
INPUTSTRPKT = 21,
MENUPKT = 6,
SYNTAXPKT = 10,
SUSPENDPKT = 17,
RESUMEPKT = 18,
BEGINDLGPKT = 19,
ENDDLGPKT = 20,
FIRSTUSERPKT = 128,
LASTUSERPKT = 255
};

enum MLTK {
MLTKERR = 0,
MLTKSTR = 34,
MLTKSYM = 35,
MLTKREAL = 42,
MLTKINT = 43,
MLTKFUNC = 70
};

typedef unsigned long long size_t;
typedef long long ssize_t;
typedef unsigned short wchar_t;
typedef long long ptrdiff_t;
typedef long long intptr_t;
typedef unsigned long long uintptr_t;
typedef signed char int8_t;
typedef signed short int int16_t;
typedef signed int int32_t;
typedef signed long long int int64_t;
typedef unsigned char uint8_t;
typedef unsigned short int uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long long int uint64_t;

typedef void *HANDLE;
struct HWND__ { int unused; }; typedef struct HWND__ *HWND;
struct HINSTANCE__ { int unused; }; typedef struct HINSTANCE__ *HINSTANCE;

typedef unsigned int wint;
typedef struct _wint{
	wint low, hi;
} wint64;
typedef long mllong32;
typedef unsigned long mlulong32;
typedef long long mlint64;
typedef unsigned long long mluint64;
typedef mlint64 mlbigint;
typedef mluint64 mlbiguint;
typedef long double    mlextended_double;

typedef int           int_ct;
typedef int_ct   mlapi_token;
typedef unsigned long      mlapi__token;
typedef mlapi__token  * mlapi__tokenp;
typedef int_ct  mlapi_packet;
typedef long mlapi_error;
typedef long mlapi__error;
typedef int_ct mlapi_result;

typedef void  * MLPointer;
typedef void  * MLENVPARAM;
typedef MLENVPARAM MLEnvironmentParameter;
typedef struct ml_environment  *MLENV;
typedef MLENV MLEnvironment;
typedef struct MLink  *MLINK;
typedef struct MLinkMark  *MLMARK;
typedef MLMARK MLINKMark;
typedef int ( *__MLProcPtr__)();
typedef void  * dev_voidp;
typedef dev_voidp dev_type;
typedef dev_type  * dev_typep;
typedef long devproc_error;
typedef unsigned long devproc_selector;
typedef struct read_buf {
	unsigned short length;
	unsigned char* ptr;
} read_buf;
typedef read_buf  * read_bufp;
typedef read_bufp  * read_bufpp;
typedef devproc_error (  * MLDeviceProcPtr) ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2);
extern devproc_error   MLDeviceMain ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2);
typedef MLDeviceProcPtr MLDeviceUPP;
typedef MLDeviceUPP dev_main_type;
typedef dev_main_type  * dev_main_typep;
typedef void* (  * MLAllocatorProcPtr) (unsigned long long );
typedef MLAllocatorProcPtr MLAllocatorUPP;
typedef void (  * MLDeallocatorProcPtr) (void*);
typedef MLDeallocatorProcPtr MLDeallocatorUPP; 
extern __MLProcPtr__  MLAllocatorCast ( MLAllocatorProcPtr f);
extern __MLProcPtr__  MLDeallocatorCast ( MLDeallocatorProcPtr f);
typedef MLAllocatorUPP MLAllocator;
typedef MLAllocator  * MLAllocatorp;
typedef MLDeallocatorUPP MLDeallocator;
typedef MLDeallocator  * MLDeallocatorp;
typedef void  * dev_world;
typedef MLINK dev_cookie;
typedef dev_world  * dev_worldp;
typedef dev_cookie  * dev_cookiep;
typedef  MLAllocatorUPP dev_allocator;
typedef  MLDeallocatorUPP dev_deallocator;
typedef dev_main_type world_main_type;
typedef unsigned long dev_mode;
typedef dev_mode  * dev_modep;
typedef unsigned long dev_options;
typedef long devyield_result;
typedef long devyield_place;
typedef long devyield_count;
typedef unsigned long devyield_sleep;
typedef struct MLYieldParams  * MLYieldParameters;
typedef struct MLYieldData{
	union {long l; double d; void  * p;} private_data[8];
}  * MLYieldDataPointer;
void MLNewYieldData ( MLYieldDataPointer ydp );
void MLFreeYieldData ( MLYieldDataPointer ydp);
MLYieldParameters MLResetYieldData ( MLYieldDataPointer ydp, devyield_place func_id);
int   MLSetYieldParameter ( MLYieldParameters yp, unsigned long selector, void* data, unsigned long* len);
int   MLYieldParameter ( MLYieldParameters yp, unsigned long selector, void* data, unsigned long* len);
devyield_sleep MLSetSleepYP ( MLYieldParameters yp, devyield_sleep sleep);
devyield_count MLSetCountYP ( MLYieldParameters yp, devyield_count count);
enum { MLSleepParameter = 1, MLCountParameter, MLPlaceParameter};
typedef int (  * MLYielderProcPtr) (MLINK mlp, MLYieldParameters yp);
typedef	MLYielderProcPtr MLDeviceYielderProcPtr;
typedef MLYielderProcPtr MLYielderUPP, MLDeviceYielderUPP;
typedef  MLYielderUPP MLYieldFunctionType;
typedef MLYielderUPP MLYieldFunctionObject;
typedef  MLYieldFunctionObject dev_yielder;
typedef dev_yielder * dev_yielderp;
typedef unsigned long dev_message;
typedef dev_message  * dev_messagep;
typedef void (  * MLHandlerProcPtr) (MLINK mlp, int m, int n);
typedef MLHandlerProcPtr MLDeviceHandlerProcPtr;
typedef MLHandlerProcPtr MLHandlerUPP, MLDeviceHandlerUPP;
typedef  MLHandlerUPP MLMessageHandlerType;
typedef MLHandlerUPP MLMessageHandlerObject;
typedef  MLMessageHandlerObject dev_msghandler;
typedef dev_msghandler * dev_msghandlerp;
extern devyield_sleep  MLSleepYP ( MLYieldParameters yp);
extern devyield_count  MLCountYP ( MLYieldParameters yp);
extern MLYieldFunctionObject  MLCreateYieldFunction ( MLEnvironment ep, MLYieldFunctionType yf, void* reserved);  
extern MLYieldFunctionType  MLDestroyYieldFunction ( MLYieldFunctionObject yfo);
extern int  MLCallYieldFunction ( MLYieldFunctionObject yfo, MLINK mlp, MLYieldParameters p);
extern MLMessageHandlerObject  MLCreateMessageHandler ( MLEnvironment ep, MLMessageHandlerType mh, void* reserved);  
extern MLMessageHandlerType  MLDestroyMessageHandler ( MLMessageHandlerObject mho);
extern void  MLCallMessageHandler ( MLMessageHandlerObject mho, MLINK mlp, int m, int n);
extern __MLProcPtr__  MLYielderCast ( MLYielderProcPtr yp);
extern __MLProcPtr__  MLHandlerCast ( MLHandlerProcPtr mh);
typedef void (  * MLSigHandlerProcPtr) (int signal);
typedef MLSigHandlerProcPtr MLSignalHandlerType;
typedef void * MLSignalHandlerObject;
typedef char  * MLParametersPointer;
typedef char MLParameters[356];
typedef void (  * MLUserProcPtr) (MLINK);
typedef MLUserProcPtr MLUserUPP;
typedef MLUserUPP MLUserFunctionType;
typedef MLUserFunctionType  * MLUserFunctionTypePointer;
typedef MLUserUPP MLUserFunction;
extern MLEnvironmentParameter  MLNewParameters (unsigned long rev, unsigned long apirev);
extern void  MLReleaseParameters (MLEnvironmentParameter ep);
extern void  MLSetAllocParameter (MLEnvironmentParameter ep, MLAllocator allocator, MLDeallocator deallocator);
extern long  MLSetThreadSafeLinksParameter (MLEnvironmentParameter ep);
extern int  MLAllocParameter (MLEnvironmentParameter ep, MLAllocator* allocator, MLDeallocator* deallocator);
extern long  MLSetResourceParameter (MLEnvironmentParameter ep, const char *path);
extern long  MLSetDeviceParameter (MLEnvironmentParameter ep, const char *devspec);
extern long  MLErrorParameter (MLEnvironmentParameter ep);
extern long  MLSetEncodingParameter (MLEnvironmentParameter ep, unsigned int etype);
extern long  MLDoNotHandleSignalParameter (MLEnvironmentParameter ep, int signum);
extern void  MLStopHandlingSignal (MLEnvironment env, int signum);
extern void  MLHandleSignal (MLEnvironment env, int signum);
extern long  MLSetEnvironmentData ( MLEnvironment env, void *cookie);
extern void *  MLEnvironmentData ( MLEnvironment env);
extern int  MLSetSignalHandler ( MLEnvironment env, int signum, void *so);
extern int  MLSetSignalHandlerFromFunction ( MLEnvironment env, int signum, MLSignalHandlerType sigfunc);
extern int  MLUnsetSignalHandler ( MLEnvironment env, int signum, MLSignalHandlerType sigfunc);
extern long  MLSetSymbolReplacement ( MLINK mlp, const char *priv, int prlen, const char *pub, int pblen);
extern int  MLClearSymbolReplacement ( MLINK mlp, long index);
extern void  MLClearAllSymbolReplacements ( MLINK mlp);
extern MLEnvironment  MLInitialize ( MLEnvironmentParameter ep);
extern void  MLDeinitialize ( MLEnvironment env);
extern void  MLVersionNumbers ( MLEnvironment env, int *inumb, int *rnumb, int *bnumb);
extern int  MLCompilerID (MLEnvironment env, const char **id);
extern void  MLReleaseCompilerID (MLEnvironment env, const char *id);
extern int  MLUCS2CompilerID (MLEnvironment env, unsigned short **id, int *length);
extern void  MLReleaseUCS2CompilerID (MLEnvironment env, unsigned short *id, int length);
extern int  MLUTF8CompilerID (MLEnvironment env, unsigned char **id, int *length);
extern void  MLReleaseUTF8CompilerID (MLEnvironment env, unsigned char *id, int length);
extern int  MLUTF16CompilerID (MLEnvironment env, unsigned short **id, int *length);
extern void  MLReleaseUTF16CompilerID (MLEnvironment env, unsigned short *id, int length);
extern int  MLUTF32CompilerID (MLEnvironment env, unsigned int **id, int *length);
extern void  MLReleaseUTF32CompilerID (MLEnvironment env, unsigned int *id, int length);
extern MLEnvironment  MLBegin (MLEnvironmentParameter ep);
extern void  MLEnd ( MLEnvironment env);
extern int  MLSetEnvIDString ( MLEnvironment ep, const char *environment_id);  
extern int  MLGetLinkedEnvIDString (MLINK mlp, const char **environment_id);  
extern void  MLReleaseEnvIDString (MLINK mlp, const char *environment_id);
extern char **  MLGetNetworkAddressList ( MLEnvironment ep, unsigned long *size );
extern void  MLReleaseNetworkAddressList ( MLEnvironment ep, char **addresses, unsigned long size);
extern char **  MLGetDomainNameList ( MLEnvironment ep, unsigned long *size );
extern void  MLReleaseDomainNameList ( MLEnvironment ep, char **dnsnames, unsigned long size);
extern int  MLGetAvailableLinkProtocolNames (MLEnvironment ep, char ***protocolNames, int *length);
extern void  MLReleaseLinkProtocolNames (MLEnvironment ep, char **protocolNames, int length);
extern int  MLGetLinksFromEnvironment (MLEnvironment ep, MLINK **links, int *length);
extern void  MLReleaseLinksFromEnvironment (MLEnvironment ep, MLINK *links, int length);
extern long  MLNumericsQuery ( MLEnvironment ep, unsigned long selector, void *p1, void *p2, long *np);
extern int  MLValid ( MLINK mlp);
extern MLINK  MLCreateLinkWithExternalProtocol ( MLEnvironment ep, dev_type dev, dev_main_type dev_main, int *errp);
extern char **  MLFilterArgv ( MLEnvironment ep, char **argv, char **argv_end);
extern long  MLFeatureString ( MLINK mlp, char *buf, long buffsize);
extern MLINK  MLOpenArgv ( MLEnvironment ep, char **argv, char **argv_end, int *errp);
extern MLINK  MLOpenArgcArgv ( MLEnvironment ep, int argc, char **argv, int *errp);
extern MLINK  MLOpenString ( MLEnvironment ep, const char *command_line, int *errp);
extern MLINK  MLLoopbackOpen ( MLEnvironment ep, int *errp);
extern int  MLStringToArgv ( const char *commandline, char *buf, char **argv, int len);
extern long  MLScanString ( char **argv, char ***argv_end, char **commandline, char **buf);
extern long  MLPrintArgv ( char *buf, char **buf_endp, char ***argvp, char **argv_end);
extern const char *  MLErrorMessage ( MLINK mlp);
extern const char *  MLErrorString ( MLEnvironment env, long err);
extern const unsigned short *  MLUCS2ErrorMessage (MLINK mlp, int *length);
extern const unsigned char *  MLUTF8ErrorMessage (MLINK mlp, int *length);
extern const unsigned short *  MLUTF16ErrorMessage (MLINK mlp, int *length);
extern const unsigned int *  MLUTF32ErrorMessage (MLINK mlp, int *length);
extern void  MLReleaseErrorMessage (MLINK mlp, const char *message);
extern void  MLReleaseUCS2ErrorMessage (MLINK mlp, const unsigned short *message, int length);
extern void  MLReleaseUTF8ErrorMessage (MLINK mlp, const unsigned char *message, int length);
extern void  MLReleaseUTF16ErrorMessage (MLINK mlp, const unsigned short *message, int length);
extern void  MLReleaseUTF32ErrorMessage (MLINK mlp, const unsigned int *message, int length);
extern MLINK  MLOpen ( int argc, char **argv);
extern MLINK  MLOpenInEnv ( MLEnvironment env, int argc, char **argv, int *errp);
extern MLINK  MLDuplicateLink ( MLINK parentmlp, const char *name, int *errp );
extern int  MLConnect ( MLINK mlp);
extern int  MLActivate ( MLINK mlp);
typedef struct feature_set* feature_setp;
extern int  MLEstablish ( MLINK mlp, feature_setp features);
extern int  MLEstablishString ( MLINK mlp, const char *features);
extern void  MLClose ( MLINK mlp);
extern void  MLSetUserData ( MLINK mlp, void* data, MLUserFunction f);
extern void*  MLUserData ( MLINK mlp, MLUserFunctionType *fp);
extern void  MLSetUserBlock ( MLINK mlp, void* userblock);
extern void*  MLUserBlock ( MLINK mlp);
extern __MLProcPtr__  MLUserCast ( MLUserProcPtr f);
extern int  MLLogStreamToFile (MLINK mlp, const char *name);
extern int  MLDisableLoggingStream (MLINK mlp);
extern int  MLEnableLoggingStream (MLINK mlp);
extern int  MLStopLoggingStreamToFile (MLINK mlp, const char *name);
extern int  MLStopLoggingStream (MLINK mlp);
extern int  MLLogFileNameForLink (MLINK mlp, const char **name);
extern void  MLReleaseLogFileNameForLink (MLINK mlp, const char *name);
extern const char *  MLName ( MLINK mlp);
extern const char *  MLLinkName ( MLINK mlp);
extern const unsigned short *  MLUCS2LinkName (MLINK mlp, int *length);
extern const unsigned char *  MLUTF8LinkName (MLINK mlp, int *length);
extern const unsigned short *  MLUTF16LinkName (MLINK mlp, int *length);
extern const unsigned int *  MLUTF32LinkName (MLINK mlp, int *length);
extern void  MLReleaseLinkName (MLINK mlp, const char *name);
extern void  MLReleaseUCS2LinkName (MLINK mlp, const unsigned short *name, int length);
extern void  MLReleaseUTF8LinkName (MLINK mlp, const unsigned char *name, int length);
extern void  MLReleaseUTF16LinkName (MLINK mlp, const unsigned short *name, int length);
extern void  MLReleaseUTF32LinkName (MLINK mlp, const unsigned int *name, int length);
extern long  MLNumber ( MLINK mlp);
extern long  MLToLinkID ( MLINK mlp);
extern MLINK  MLFromLinkID ( MLEnvironment ep, long n);
extern char *  MLSetName ( MLINK mlp, const char *name);
extern void*  MLInit ( MLAllocator alloc, MLDeallocator dealloc, void* enclosing_environment);
extern void  MLDeinit ( void* env);
extern void*  MLEnclosingEnvironment ( void* ep);
extern void*  MLinkEnvironment ( MLINK mlp);
extern void  MLEnableLinkLock (MLINK mlp);
extern void  MLDisableLinkLock (MLINK mlp);
extern MLEnvironment  MLLinkEnvironment (MLINK mlp);
extern int  MLIsLinkLoopback (MLINK mlp);
extern MLYieldFunctionObject  MLDefaultYieldFunction ( MLEnvironment env);
extern int  MLSetDefaultYieldFunction ( MLEnvironment env, MLYieldFunctionObject yf);
typedef void * MLLinkServer;
typedef void (*MLNewLinkCallbackFunction)(MLLinkServer server, MLINK link);
extern MLLinkServer  MLNewLinkServer (MLEnvironment env, void *context, int *error);
extern MLLinkServer  MLNewLinkServerWithPort (MLEnvironment env, unsigned short port, void *context, int *error);
extern MLLinkServer  MLNewLinkServerWithPortAndInterface (MLEnvironment env, unsigned short port, const char *iface, void *context, int *error);
extern void  MLShutdownLinkServer (MLLinkServer server);
extern void  MLRegisterCallbackFunctionWithLinkServer (MLLinkServer server, MLNewLinkCallbackFunction function);
extern MLINK  MLWaitForNewLinkFromLinkServer (MLLinkServer server, int *error);
extern unsigned short  MLPortFromLinkServer (MLLinkServer server, int *error);
extern const char *  MLInterfaceFromLinkServer (MLLinkServer server, int *error);
extern void *  MLContextFromLinkServer (MLLinkServer server, int *error);
extern void  MLReleaseInterfaceFromLinkServer (MLLinkServer server, const char *iface);
typedef void * MLServiceRef;
typedef void (*MLBrowseCallbackFunction)(MLEnvironment env, MLServiceRef ref, int flag,
	const char *serviceName, void *context);
extern int  MLBrowseForLinkServices (MLEnvironment env, MLBrowseCallbackFunction callbackFunction, const char *serviceProtocol, const char *domain, void *context, MLServiceRef *ref);
extern void  MLStopBrowsingForLinkServices (MLEnvironment env, MLServiceRef ref);
typedef void (*MLResolveCallbackFunction)(MLEnvironment env, MLServiceRef ref, const char *serviceName,
	const char *linkName, const char *protocol, int options, void *context);
extern int  MLResolveLinkService (MLEnvironment env, MLResolveCallbackFunction, const char *serviceProtocol, const char *serviceName, void *context, MLServiceRef *ref);
extern void  MLStopResolvingLinkService (MLEnvironment env, MLServiceRef ref);
typedef void (*MLRegisterCallbackFunction)(MLEnvironment env, MLServiceRef ref, int flag, const char *serviceName,
	void *context);
extern MLINK  MLRegisterLinkServiceWithPortAndHostname (MLEnvironment env, const char *serviceProtocol, const char *serviceName, unsigned short port, const char *hostname, MLRegisterCallbackFunction function, const char *domain, void *context, MLServiceRef *ref, int *error);
extern MLINK  MLRegisterLinkServiceWithHostname (MLEnvironment env, const char *serviceProtocol, const char *serviceName, const char *hostname, MLRegisterCallbackFunction function, const char *domain, void *context, MLServiceRef *ref, int *error);
extern MLINK  MLRegisterLinkService (MLEnvironment env, const char *serviceProtocol, const char *serviceName, MLRegisterCallbackFunction function, const char *domain, void *context, MLServiceRef *, int *error);
extern MLINK  MLRegisterLinkServiceUsingLinkProtocol (MLEnvironment env, const char *serviceProtocol, const char *serviceName, unsigned short port, const char *hostname, const char *protocol, MLRegisterCallbackFunction function, const char *domain, void *context, MLServiceRef *ref, int *error);
extern void  MLRegisterLinkServiceFromLinkServer (MLEnvironment env, const char *serviceProtocol, const char *serviceName, MLLinkServer server, MLRegisterCallbackFunction function, const char *domain, void *context, MLServiceRef *ref, int *error);
extern void  MLStopRegisteringLinkService (MLEnvironment env, MLServiceRef ref);
extern void  MLStopRegisteringLinkServiceForLink (MLEnvironment env, MLINK link, MLServiceRef ref);
extern const char *  MLServiceProtocolFromReference (MLEnvironment env, MLServiceRef ref);
extern int  MLError ( MLINK mlp);
extern int  MLClearError ( MLINK mlp);
extern int  MLSetError ( MLINK mlp, int err);
enum {	MLTerminateMessage = 1, MLInterruptMessage, MLAbortMessage,
	MLEndPacketMessage, MLSynchronizeMessage, MLImDyingMessage,
	MLWaitingAcknowledgment, MLMarkTopLevelMessage, MLLinkClosingMessage,
	MLAuthenticateFailure, MLSuspendActivitiesMessage, MLResumeActivitiesMessage,
	MLFirstUserMessage = 128, MLLastUserMessage = 255 };
typedef unsigned long devinfo_selector;
extern int  MLPutMessage ( MLINK mlp, int msg);
extern int  MLGetMessage ( MLINK mlp, int *mp, int *np);
extern int  MLMessageReady ( MLINK mlp);
extern int  MLPutMessageWithArg ( MLINK mlp, int msg, int arg);
extern MLMessageHandlerObject  MLGetMessageHandler ( MLINK mlp);
extern MLMessageHandlerObject  MLMessageHandler ( MLINK mlp);
extern MLYieldFunctionObject  MLGetYieldFunction ( MLINK mlp);
extern MLYieldFunctionObject  MLYieldFunction ( MLINK mlp);
extern int  MLSetMessageHandler ( MLINK mlp, MLMessageHandlerObject h);
extern int  MLSetYieldFunction ( MLINK mlp, MLYieldFunctionObject yf);
extern int  MLDeviceInformation ( MLINK mlp, devinfo_selector selector, void* buf, long *buflen);
extern int  MLLowLevelDeviceName (MLINK mlp, const char **name);
extern void  MLReleaseLowLevelDeviceName (MLINK mlp, const char *name);
extern int  MLGetNext ( MLINK mlp);
extern int  MLGetNextRaw ( MLINK mlp);
extern int  MLGetType ( MLINK mlp);
extern int  MLGetRawType ( MLINK mlp);
extern int  MLGetRawData ( MLINK mlp, unsigned char *data, int size, int *gotp);
extern int  MLGetData ( MLINK mlp, char *data, int size, int *gotp);
extern int  MLGetArgCount ( MLINK mlp, int *countp);
extern int  MLGetRawArgCount ( MLINK mlp, int *countp);
extern int  MLBytesToGet ( MLINK mlp, int *leftp);
extern int  MLRawBytesToGet ( MLINK mlp, int *leftp);
extern int  MLExpressionsToGet ( MLINK mlp, int *countp);
extern int  MLNewPacket ( MLINK mlp);
extern int  MLTakeLast ( MLINK mlp, int eleft);
extern int  MLFill ( MLINK mlp);
extern int  MLPutNext ( MLINK mlp, int tok);
extern int  MLPutType ( MLINK mlp, int tok);
extern int  MLPutRawSize ( MLINK mlp, int size);
extern int  MLPutRawData ( MLINK mlp, const unsigned char *data, int len);
extern int  MLPutArgCount ( MLINK mlp, int argc);
extern int  MLPutComposite ( MLINK mlp, int argc);
extern int  MLBytesToPut ( MLINK mlp, int *leftp);
extern int  MLEndPacket ( MLINK mlp);
extern int  MLFlush ( MLINK mlp);
typedef unsigned long decoder_mask;
extern int  MLGetBinaryNumber ( MLINK mlp, void *np, long type);
extern int  MLGetShortInteger ( MLINK mlp, short *hp);
extern int  MLGetInteger ( MLINK mlp, int *ip);
extern int  MLGetLongInteger ( MLINK mlp, long *lp);
extern int  MLGetInteger16 ( MLINK mlp, short *hp);
extern int  MLGetInteger32 ( MLINK mlp, int *ip);
extern int  MLGetInteger64 ( MLINK mlp, mlint64 *wp);
extern int  MLGetInteger8 (MLINK mlp, unsigned char *cp);
extern int  MLGetFloat ( MLINK mlp, float *fp);
extern int  MLGetDouble ( MLINK mlp, double *dp);
extern int  MLGetReal ( MLINK mlp, double *dp);
extern int  MLGetLongDouble ( MLINK mlp, mlextended_double *xp);
extern int  MLGetReal32 ( MLINK mlp, float *fp);
extern int  MLGetReal64 ( MLINK mlp, double *dp);
extern int  MLGetReal128 ( MLINK mlp, mlextended_double *dp);
extern int  MLGet8BitCharacters ( MLINK mlp, long *chars_left, unsigned char *buf, long cardof_buf, long *got, long missing);
extern int  MLGet7BitCharacters ( MLINK mlp, long *chars_left, char *buf, long cardof_buf, long *got);
extern int  MLGetUCS2Characters ( MLINK mlp, int *chars_left, unsigned short *buf, int cardof_buf, int *got);
extern int  MLGetUTF8Characters ( MLINK mlp, int *chars_left, unsigned char *buf, int cardof_buf, int *got);
extern int  MLGetUTF16Characters ( MLINK mlp, int *chars_left, unsigned short *buf, int cardof_buf, int *got);
extern int  MLGetUTF32Characters ( MLINK mlp, int *chars_left, unsigned int *buf, int cardof_buf, int *got);
extern int  MLGetByteString ( MLINK mlp, const unsigned char **sp, int *lenp, long missing);
extern int  MLGetString ( MLINK mlp, const char **sp);
extern int  MLGetUCS2String ( MLINK mlp, const unsigned short **sp, int *lenp);
extern int  MLGetUTF8String ( MLINK mlp, const unsigned char **sp, int *bytes, int *chars);
extern int  MLGetUTF16String ( MLINK mlp, const unsigned short **sp, int *ncodes, int *chars);
extern int  MLGetUTF32String ( MLINK mlp, const unsigned int **sp, int *len);
extern int  MLGetNumberAsByteString ( MLINK mlp, const unsigned char **sp, long *lenp, long missing);
extern int  MLGetNumberAsString ( MLINK mlp, const char **sp);
extern int  MLGetNumberAsUCS2String ( MLINK mlp, const unsigned short **sp, int *lenp);
extern int  MLGetNumberAsUTF8String ( MLINK mlp, const unsigned char **sp, int *bytes, int *chars);
extern int  MLGetNumberAsUTF16String ( MLINK mlp, const unsigned short **sp, int *ncodes, int *chars);
extern int  MLGetNumberAsUTF32String ( MLINK mlp, const unsigned int **sp, int *lenp);
extern void  MLReleaseUCS2String ( MLINK mlp, const unsigned short *s, int len);
extern void  MLReleaseUTF8String ( MLINK mlp, const unsigned char *s, int len);
extern void  MLReleaseUTF16String ( MLINK mlp, const unsigned short *s, int len);
extern void  MLReleaseUTF32String ( MLINK mlp, const unsigned int *s, int len);
extern void  MLReleaseByteString ( MLINK mlp, const unsigned char * s, int len);
extern void  MLReleaseString ( MLINK mlp, const char *s);
extern int  MLTestString ( MLINK mlp, const char *name);
extern int  MLTestUCS2String ( MLINK mlp, const unsigned short *name, int length);
extern int  MLTestUTF8String ( MLINK mlp, const unsigned char *name, int length);
extern int  MLTestUTF16String ( MLINK mlp, const unsigned short *name, int length);
extern int  MLTestUTF32String ( MLINK mlp, const unsigned int *name, int length);
extern int  MLGetByteSymbol ( MLINK mlp, const unsigned char ** sp, int *lenp, long missing);
extern int  MLGetSymbol ( MLINK mlp, const char ** sp);
extern int  MLGetUCS2Symbol ( MLINK mlp, const unsigned short **sp, int *lenp);
extern int  MLGetUTF8Symbol ( MLINK mlp, const unsigned char **sp, int *bytes, int *chars);
extern int  MLGetUTF16Symbol ( MLINK mlp, const unsigned short **sp, int *ncodes, int *chars);
extern int  MLGetUTF32Symbol ( MLINK mlp, const unsigned int **sp, int *lenp);
extern void  MLReleaseUCS2Symbol ( MLINK mlp, const unsigned short *s, int len);
extern void  MLReleaseUTF8Symbol ( MLINK mlp, const unsigned char *s, int len);
extern void  MLReleaseUTF16Symbol ( MLINK mlp, const unsigned short *s, int len);
extern void  MLReleaseUTF32Symbol ( MLINK mlp, const unsigned int *s, int len);
extern void  MLReleaseByteSymbol ( MLINK mlp, const unsigned char * s, int len);
extern void  MLReleaseSymbol ( MLINK mlp, const char *s);
extern int  MLTestSymbol ( MLINK mlp, const char *name);
extern int  MLTestUCS2Symbol ( MLINK mlp, const unsigned short *name, int length);
extern int  MLTestUTF8Symbol ( MLINK mlp, const unsigned char *name, int length);
extern int  MLTestUTF16Symbol ( MLINK mlp, const unsigned short *name, int length);
extern int  MLTestUTF32Symbol ( MLINK mlp, const unsigned int *name, int length);
extern int  MLGetFunction ( MLINK mlp, const char **sp, int *countp);
extern int  MLGetUCS2Function ( MLINK mlp, const unsigned short **sp, int *length, int *countp);
extern int  MLGetUTF8Function ( MLINK mlp, const unsigned char **sp, int *length, int *countp);
extern int  MLGetUTF16Function ( MLINK mlp, const unsigned short **sp, int *length, int *countp);
extern int  MLGetUTF32Function ( MLINK mlp, const unsigned int **sp, int *length, int *countp);
extern int  MLCheckFunction ( MLINK mlp, const char *s, long *countp);
extern int  MLCheckFunctionWithArgCount ( MLINK mlp, const char *s, long *countp);
extern int  MLTestHead ( MLINK mlp, const char *s, int *countp);
extern int  MLTestHeadWithArgCount (MLINK mlp, const char *s, int *countp);
extern int  MLTestUCS2HeadWithArgCount (MLINK mlp, const unsigned short *s, int length, int *countp);
extern int  MLTestUTF16HeadWithArgCount (MLINK mlp, const unsigned short *s, int length, int *countp);
extern int  MLTestUTF32HeadWithArgCount (MLINK mlp, const unsigned int *s, int length, int *countp);
extern int  MLTestUTF8HeadWithArgCount (MLINK mlp, const unsigned char *s, int length, int *countp);
extern int  MLTestUCS2Head ( MLINK mlp, const unsigned short *s, int length, int *countp);
extern int  MLTestUTF8Head ( MLINK mlp, const unsigned char *s, int length, int *countp);
extern int  MLTestUTF16Head ( MLINK mlp, const unsigned short *s, int length, int *countp);
extern int  MLTestUTF32Head ( MLINK mlp, const unsigned int *s, int length, int *countp);
extern int  MLPutBinaryNumber ( MLINK mlp, void *np, long type);
extern int  MLPutShortInteger ( MLINK mlp, int h);
extern int  MLPutInteger ( MLINK mlp, int i);
extern int  MLPutLongInteger ( MLINK mlp, long l);
extern int  MLPutInteger16 ( MLINK mlp, int h);
extern int  MLPutInteger32 ( MLINK mlp, int i);
extern int  MLPutInteger64 ( MLINK mlp, mlint64 w);
extern int  MLPutInteger8 (MLINK mlp, unsigned char i);
extern int  MLPutFloat ( MLINK mlp, double f);
extern int  MLPutDouble ( MLINK mlp, double d);
extern int  MLPutReal ( MLINK mlp, double d);
extern int  MLPutLongDouble ( MLINK mlp, mlextended_double x);
extern int  MLPutReal32 ( MLINK mlp, double f);
extern int  MLPutReal64 ( MLINK mlp, double d);
extern int  MLPutReal128 ( MLINK mlp, mlextended_double x);
extern int  MLPut8BitCharacters ( MLINK mlp, long chars_left, const unsigned char *bytes, long nbytes);
extern int  MLPut7BitCount ( MLINK mlp, long count, long size);
extern int  MLPut7BitCharacters ( MLINK mlp, long chars_left, const char *bytes, long nbytes, long nchars_now);
extern int  MLPutUCS2Characters ( MLINK mlp, int chars_left, const unsigned short *codes, int ncodes);
extern int  MLPutUTF8Characters ( MLINK mlp, int chars_left, const unsigned char *codes, int ncodes);
extern int  MLPutUTF16Characters ( MLINK mlp, int chars_left, const unsigned short *codes, int ncodes);
extern int  MLPutUTF32Characters ( MLINK mlp, int chars_left, const unsigned int *codes, int ncodes);
extern int  MLPutByteString ( MLINK mlp, const unsigned char *s, long len);
extern int  MLPutString ( MLINK mlp, const char *s);
extern int  MLPutUCS2String ( MLINK mlp, const unsigned short *s, int len);
extern int  MLPutUTF8String ( MLINK mlp, const unsigned char *s, int len);
extern int  MLPutUTF16String ( MLINK mlp, const unsigned short *s, int len);
extern int  MLPutUTF32String ( MLINK mlp, const unsigned int *s, int len);
extern int  MLPutRealNumberAsString ( MLINK mlp, const char *s);
extern int  MLPutRealNumberAsByteString ( MLINK mlp, const unsigned char *s);
extern int  MLPutRealNumberAsUCS2String ( MLINK mlp, const unsigned short *s);
extern int  MLPutRealNumberAsUTF8String ( MLINK mlp, const unsigned char *s, int nbytes);
extern int  MLPutRealNumberAsUTF16String ( MLINK mlp, const unsigned short *s, int ncodes);
extern int  MLPutRealNumberAsUTF32String ( MLINK mlp, const unsigned int *s, int nchars);
extern int  MLPutSize ( MLINK mlp, int size);
extern int  MLPutData ( MLINK mlp, const char *buff, int len);
extern int  MLPutByteSymbol ( MLINK mlp, const unsigned char *s, long len);
extern int  MLPutSymbol ( MLINK mlp, const char *s);
extern int  MLPutUCS2Symbol ( MLINK mlp, const unsigned short *s, int len);
extern int  MLPutUTF8Symbol ( MLINK mlp, const unsigned char *s, int len);
extern int  MLPutUTF16Symbol ( MLINK mlp, const unsigned short *s, int len);
extern int  MLPutUTF32Symbol ( MLINK mlp, const unsigned int *s, int len);
extern int  MLPutFunction ( MLINK mlp, const char *s, int argc);
extern int  MLPutUCS2Function ( MLINK mlp, const unsigned short *s, int length, int argn);
extern int  MLPutUTF8Function ( MLINK mlp, const unsigned char *s, int length, int argn);
extern int  MLPutUTF16Function ( MLINK mlp, const unsigned short *s, int length, int argn);
extern int  MLPutUTF32Function ( MLINK mlp, const unsigned int *s, int length, int argn);
typedef struct {
	const char *str;
	const char *end;
} MLStringPosition;
typedef MLStringPosition  * MLStringPositionPointer;
typedef struct {
	unsigned char *cc;
	int  mode;
	int  more;
	unsigned char *head;
} MLOldStringPosition;
typedef MLOldStringPosition  * MLOldStringPositionPointer;
extern long  MLCharacterOffset ( const char **startp, const char *end, long n);
extern long  MLStringCharacter ( const char * start, const char *end);
extern long  MLNextCharacter ( const char **startp, const char *end);
extern long  MLNextCharacterFromStringWithLength (const char *str, long *indexp, long len);
extern long  MLConvertNewLine ( char **sp);
extern long  MLConvertCharacter ( unsigned long ch, char **sp);
extern long  MLConvertByteString ( unsigned char *codes, long len, char **strp, char *str_end);
extern long  MLConvertByteStringNL ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl);
extern long  MLConvertDoubleByteString ( unsigned char *codes, long len, char **strp, char *str_end);
extern long  MLConvertDoubleByteStringNL ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl);
extern long  MLConvertUCS2String ( unsigned short *codes, long len, char **strp, char *str_end);
extern long  MLConvertUCS2StringNL ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl);
extern long  MLConvertUTF8String ( unsigned char *codes, long len, char **strp, char *str_end);
extern long  MLConvertUTF8StringNL ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl);
extern long  MLConvertUTF16String ( unsigned short *codes, long len, char **strp, char *str_end);
extern long  MLConvertUTF16StringNL ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl);
extern long  MLConvertUTF32String ( unsigned int *codes, long len, char **strp, char *str_end);
extern long  MLConvertUTF32StringNL ( unsigned int *codes, long len, char **strp, char *str_end, unsigned long nl);
extern const char *  MLStringFirstPosFun ( const char *s, MLStringPositionPointer p);
extern int  MLOldPutCharToString ( unsigned int ch, char **sp);
extern unsigned char *  MLOldStringNextPosFun ( MLOldStringPositionPointer p);
extern unsigned char *  MLOldStringFirstPosFun ( char *s, MLOldStringPositionPointer p);
extern unsigned int  MLOldStringCharFun ( MLOldStringPositionPointer p);
extern long  MLOldConvertByteString ( unsigned char *codes, long len, char **strp, char *str_end);
extern long  MLOldConvertUCS2String ( unsigned short *codes, long len, char **strp, char *str_end);
extern long  MLCharOffset ( const char **startp, const char *end, long n, int more);
extern long  MLNextChar ( const char **startp, const char *end, int more, int useSurrogates, int *wasSurrogatePair);
typedef struct array_meter  * array_meterp;
typedef array_meterp  * array_meterpp;
extern int  MLPutArray ( MLINK mlp, array_meterp meterp);
extern int  MLPutBinaryNumberArrayData ( MLINK mlp, array_meterp meterp, const void * datap, long count, long type);
extern int  MLPutByteArrayData ( MLINK mlp, array_meterp meterp, const unsigned char *datap, long count);
extern int  MLPutShortIntegerArrayData ( MLINK mlp, array_meterp meterp, const short * datap, long count);
extern int  MLPutIntegerArrayData ( MLINK mlp, array_meterp meterp, const int * datap, long count);
extern int  MLPutLongIntegerArrayData ( MLINK mlp, array_meterp meterp, const long * datap, long count);
extern int  MLPutInteger8ArrayData ( MLINK mlp, array_meterp meterp, const unsigned char * datap, int count);
extern int  MLPutInteger16ArrayData ( MLINK mlp, array_meterp meterp, const short * datap, int count);
extern int  MLPutInteger32ArrayData ( MLINK mlp, array_meterp meterp, const int * datap, int count);
extern int  MLPutInteger64ArrayData ( MLINK mlp, array_meterp meterp, const mlint64 * datap, int count);
extern int  MLPutFloatArrayData ( MLINK mlp, array_meterp meterp, const float * datap, long count);
extern int  MLPutDoubleArrayData ( MLINK mlp, array_meterp meterp, const double *datap, long count);
extern int  MLPutLongDoubleArrayData ( MLINK mlp, array_meterp meterp, const mlextended_double *datap, long count);
extern int  MLPutReal32ArrayData ( MLINK mlp, array_meterp meterp, const float * datap, int count);
extern int  MLPutReal64ArrayData ( MLINK mlp, array_meterp meterp, const double *datap, int count);
extern int  MLPutReal128ArrayData ( MLINK mlp, array_meterp meterp, const mlextended_double *datap, int count);
extern int  MLPutBinaryNumberArray ( MLINK mlp, const void * data, const long *dimp, const char **heads, long depth, long type);
extern int  MLPutByteArray ( MLINK mlp, const unsigned char *data, const int *dims, const char **heads, int depth);
extern int  MLPutShortIntegerArray ( MLINK mlp, const short * data, const long *dims, const char **heads, long depth);
extern int  MLPutIntegerArray ( MLINK mlp, const int * data, const long *dims, const char **heads, long depth);
extern int  MLPutLongIntegerArray ( MLINK mlp, const long * data, const long *dims, const char **heads, long depth);
extern int  MLPutInteger8Array ( MLINK mlp, const unsigned char *data, const int *dims, const char **heads, int depth);
extern int  MLPutInteger16Array ( MLINK mlp, const short * data, const int *dims, const char **heads, int depth);
extern int  MLPutInteger32Array ( MLINK mlp, const int * data, const int *dims, const char **heads, int depth);
extern int  MLPutInteger64Array ( MLINK mlp, const mlint64 * data, const int *dims, const char **heads, int depth);
extern int  MLPutFloatArray ( MLINK mlp, const float * data, const long *dims, const char **heads, long depth);
extern int  MLPutDoubleArray ( MLINK mlp, const double *data, const long *dims, const char **heads, long depth);
extern int  MLPutRealArray ( MLINK mlp, const double *data, const long *dims, const char **heads, long depth);
extern int  MLPutLongDoubleArray ( MLINK mlp, const mlextended_double *data, const long *dims, const char **heads, long depth);
extern int  MLPutReal32Array ( MLINK mlp, const float * data, const int *dims, const char **heads, int depth);
extern int  MLPutReal64Array ( MLINK mlp, const double *data, const int *dims, const char **heads, int depth);
extern int  MLPutReal128Array ( MLINK mlp, const mlextended_double *data, const int *dims, const char **heads, int depth);
extern int  MLPutBinaryNumberList ( MLINK mlp, const void * data, long count, long type);
extern int  MLPutIntegerList ( MLINK mlp, const int * data, long count);
extern int  MLPutRealList ( MLINK mlp, const double *data, long count);
extern int  MLPutInteger8List ( MLINK mlp, const unsigned char *data, int count);
extern int  MLPutInteger16List ( MLINK mlp, const short * data, int count);
extern int  MLPutInteger32List ( MLINK mlp, const int * data, int count);
extern int  MLPutInteger64List ( MLINK mlp, const mlint64 * data, int count);
extern int  MLPutReal32List ( MLINK mlp, const float * data, int count);
extern int  MLPutReal64List ( MLINK mlp, const double *data, int count);
extern int  MLPutReal128List ( MLINK mlp, const mlextended_double *data, int count);
extern int  MLPutArrayType ( MLINK mlp, MLINK heads, long depth, array_meterpp meterpp);
extern int  MLReleasePutArrayState ( MLINK mlp, MLINK heads, array_meterp meterp);
extern int  MLPutArrayLeaves ( MLINK mlp, MLINK heads, array_meterp meterp, MLINK leaves, long count);
extern int  MLPutBinaryNumberArrayDataWithHeads ( MLINK mlp, MLINK heads, array_meterp meterp, const void *datap, long count, long type);
extern int  MLGetArrayDimensions ( MLINK mlp, array_meterp meterp);
extern int  MLGetArrayType ( MLINK mlp, array_meterp meterp);
extern int  MLGetBinaryNumberList ( MLINK mlp, void **datap, long *countp, long type);
extern int  MLGetIntegerList ( MLINK mlp, int **datap, long *countp);
extern int  MLGetRealList ( MLINK mlp, double **datap, long *countp);
extern int  MLGetInteger16List ( MLINK mlp, short ** datap, int *countp);
extern int  MLGetInteger32List ( MLINK mlp, int ** datap, int *countp);
extern int  MLGetInteger64List ( MLINK mlp, mlint64 ** datap, int *countp);
extern int  MLGetReal32List ( MLINK mlp, float ** datap, int *countp);
extern int  MLGetReal64List ( MLINK mlp, double ** datap, int *countp);
extern int  MLGetReal128List ( MLINK mlp, mlextended_double ** datap, int *countp);
extern void  MLReleaseIntegerList ( MLINK mlp, int *data, long count);
extern void  MLReleaseRealList ( MLINK mlp, double *data, long count);
extern void  MLReleaseBinaryNumberList ( MLINK mlp, void *data, int count, long type);
extern void  MLReleaseInteger16List ( MLINK mlp, short *data, int count);
extern void  MLReleaseInteger32List ( MLINK mlp, int *data, int count);
extern void  MLReleaseInteger64List ( MLINK mlp, mlint64 *data, int count);
extern void  MLReleaseReal32List ( MLINK mlp, float *data, int count);
extern void  MLReleaseReal64List ( MLINK mlp, double *data, int count);
extern void  MLReleaseReal128List ( MLINK mlp, mlextended_double *data, int count);
extern int  MLGetBinaryNumberArrayData ( MLINK mlp, array_meterp meterp, void *datap, long count, long type);
extern int  MLGetByteArrayData ( MLINK mlp, array_meterp meterp, unsigned char * datap, long count);
 
extern int  MLGetShortIntegerArrayData ( MLINK mlp, array_meterp meterp, short * datap, long count);
extern int  MLGetIntegerArrayData ( MLINK mlp, array_meterp meterp, int * datap, long count);
extern int  MLGetLongIntegerArrayData ( MLINK mlp, array_meterp meterp, long * datap, long count);
extern int  MLGetInteger16ArrayData ( MLINK mlp, array_meterp meterp, short * datap, int count);
extern int  MLGetInteger32ArrayData ( MLINK mlp, array_meterp meterp, int * datap, int count);
extern int  MLGetInteger64ArrayData ( MLINK mlp, array_meterp meterp, mlint64 * datap, int count);
extern int  MLGetFloatArrayData ( MLINK mlp, array_meterp meterp, float *datap, long count);
extern int  MLGetDoubleArrayData ( MLINK mlp, array_meterp meterp, double *datap, long count);
extern int  MLGetLongDoubleArrayData ( MLINK mlp, array_meterp meterp, mlextended_double *datap, long count);
extern int  MLGetReal32ArrayData ( MLINK mlp, array_meterp meterp, float *datap, int count);
extern int  MLGetReal64ArrayData ( MLINK mlp, array_meterp meterp, double *datap, int count);
extern int  MLGetReal128ArrayData ( MLINK mlp, array_meterp meterp, mlextended_double *datap, int count);
extern int  MLGetInteger8List (MLINK mlp, unsigned char **datap, int *countp);
extern int  MLGetInteger8ArrayData (MLINK mlp, array_meterp meterp, unsigned char *datap, int count);
extern void  MLReleaseInteger8List (MLINK mlp, unsigned char *data, int count);
extern int  MLGetArrayTypeWithDepthAndLeafType ( MLINK mlp, MLINK heads, array_meterpp meterpp, long *depthp, mlapi__token *leaf_tokp);
extern int  MLGetBinaryNumberArrayDataWithHeads ( MLINK mlp, MLINK heads, array_meterp meterp, void *datap, long *countp, long type);
extern void  MLReleaseGetArrayState ( MLINK mlp, MLINK heads, array_meterp meterp);
extern int  MLGetBinaryNumberArrayWithLeafType ( MLINK mlp, void **datap, long **dimpp, char ***headsp, long *depthp, long type, mlapi__token *leaf_tokp);
extern int  MLGetBinaryNumberArray ( MLINK mlp, void ** datap, long **dimpp, char ***headsp, long *depthp, long type);
extern int  MLGetByteArray ( MLINK mlp, unsigned char ** datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetShortIntegerArray ( MLINK mlp, short ** datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetIntegerArray ( MLINK mlp, int ** datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetLongIntegerArray ( MLINK mlp, long ** datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetInteger16Array ( MLINK mlp, short ** datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetInteger32Array ( MLINK mlp, int ** datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetInteger64Array ( MLINK mlp, mlint64 ** datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetInteger8Array (MLINK mlp, unsigned char **datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetFloatArray ( MLINK mlp, float ** datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetDoubleArray ( MLINK mlp, double **datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetRealArray ( MLINK mlp, double **datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetLongDoubleArray ( MLINK mlp, mlextended_double **datap, long **dimsp, char ***headsp, long *depthp);
extern int  MLGetReal32Array ( MLINK mlp, float ** datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetReal64Array ( MLINK mlp, double **datap, int **dimsp, char ***headsp, int *depthp);
extern int  MLGetReal128Array ( MLINK mlp, mlextended_double **datap, int **dimsp, char ***headsp, int *depthp);
extern void  MLReleaseShortIntegerArray ( MLINK mlp, short * data, long *dims, char **heads, long depth);
extern void  MLReleaseIntegerArray ( MLINK mlp, int * data, long *dims, char **heads, long depth);
extern void  MLReleaseLongIntegerArray ( MLINK mlp, long * data, long *dims, char **heads, long depth);
extern void  MLReleaseBinaryNumberArray ( MLINK mlp, void * data, int *dimp, char **heads, int len, long type);
extern void  MLReleaseByteArray ( MLINK mlp, unsigned char *data, int *dims, char **heads, int depth);
extern void  MLReleaseInteger16Array ( MLINK mlp, short * data, int *dims, char **heads, int depth);
extern void  MLReleaseInteger32Array ( MLINK mlp, int * data, int *dims, char **heads, int depth);
extern void  MLReleaseInteger64Array ( MLINK mlp, mlint64 * data, int *dims, char **heads, int depth);
extern void  MLReleaseInteger8Array (MLINK mlp, unsigned char *data, int *dimp, char **heads, int depth);
extern void  MLReleaseFloatArray ( MLINK mlp, float * data, long *dims, char **heads, long depth);
extern void  MLReleaseDoubleArray ( MLINK mlp, double *data, long *dims, char **heads, long depth);
extern void  MLReleaseRealArray ( MLINK mlp, double *data, long *dims, char **heads, long depth);
extern void  MLReleaseReal32Array ( MLINK mlp, float * data, int *dims, char **heads, int depth);
extern void  MLReleaseReal64Array ( MLINK mlp, double *data, int *dims, char **heads, int depth);
extern void  MLReleaseReal128Array ( MLINK mlp, mlextended_double *data, int *dims, char **heads, int depth);
extern void  MLReleaseLongDoubleArray ( MLINK mlp, mlextended_double *data, long *dims, char **heads, long depth);
enum MLUnicodeContainerType
{
	UCS2ContainerType,
	UTF8ContainerType,
	UTF16ContainerType,
	UTF32ContainerType
};
typedef struct _MLUnicodeContainer
{
	union _pointer
	{
		unsigned short *ucs2;
		unsigned char *utf8;
		unsigned short *utf16;
		unsigned int *utf32;
	} pointer;
	int length;
	enum MLUnicodeContainerType type;
} MLUnicodeContainer;
extern MLUnicodeContainer *  MLNewUnicodeContainer (const void *string, int length, enum MLUnicodeContainerType type);
extern void  MLReleaseUnicodeContainer (MLUnicodeContainer *string);
extern MLINKMark  MLCreateMark ( MLINK mlp);
extern MLINKMark  MLSeekToMark ( MLINK mlp, MLINKMark mark, int index);
extern MLINKMark  MLSeekMark ( MLINK mlp, MLINKMark mark, int index);
extern void  MLDestroyMark ( MLINK mlp, MLINKMark mark);
extern int  MLTransferExpression ( MLINK dmlp, MLINK smlp);
extern int  MLTransferToEndOfLoopbackLink ( MLINK dmlp, MLINK smlp);
extern int  MLForwardReset ( MLINK mlp, unsigned long marker);
extern int  MLAlign ( MLINK lmlp, MLINK rmlp);
extern int  MLNextPacket ( MLINK mlp);
typedef long long mldlg_result;
typedef mldlg_result (  * MLAlertProcPtr) ( MLEnvironment env, const char *message);
typedef mldlg_result (  * MLRequestProcPtr) ( MLEnvironment env, const char *prompt, char *response, long sizeof_response);
typedef mldlg_result (  * MLConfirmProcPtr) ( MLEnvironment env, const char *question, mldlg_result default_answer);
typedef mldlg_result (  * MLRequestArgvProcPtr) ( MLEnvironment env, char **argv, long cardof_argv, char *buf, long sizeof_buf);
typedef mldlg_result (  * MLRequestToInteractProcPtr) ( MLEnvironment env, mldlg_result wait_for_permission);
typedef mldlg_result (  * MLDialogProcPtr) ( MLEnvironment env);
typedef MLDialogProcPtr MLDialogUPP;
typedef MLAlertProcPtr MLAlertUPP;
typedef MLRequestProcPtr MLRequestUPP;
typedef MLConfirmProcPtr MLConfirmUPP;
typedef MLRequestArgvProcPtr MLRequestArgvUPP;
typedef MLRequestToInteractProcPtr MLRequestToInteractUPP;
typedef MLAlertUPP MLAlertFunctionType;
typedef MLRequestUPP MLRequestFunctionType;
typedef MLConfirmUPP MLConfirmFunctionType;
typedef MLRequestArgvUPP MLRequestArgvFunctionType;
typedef MLRequestToInteractUPP MLRequestToInteractFunctionType;
typedef MLDialogUPP MLDialogFunctionType;
enum {	MLAlertFunction = 1, MLRequestFunction, MLConfirmFunction,
	MLRequestArgvFunction, MLRequestToInteractFunction };
extern mldlg_result   MLAlert_win ( MLEnvironment ep, const char *alertstr);
extern mldlg_result   MLRequest_win ( MLEnvironment ep, const char *prompt, char *response, long n);
extern mldlg_result   MLConfirm_win ( MLEnvironment ep, const char *okcancelquest, mldlg_result default_answer);
extern mldlg_result   MLPermit_win ( MLEnvironment ep, mldlg_result wait);
extern mldlg_result   default_request_argv ( MLEnvironment ep, char **argv, long len, char *buff, long size);
extern mldlg_result  MLAlert ( MLEnvironment env, const char *message);
extern mldlg_result  MLRequest ( MLEnvironment env, const char *prompt, char *response, long sizeof_response);  
extern mldlg_result  MLConfirm ( MLEnvironment env, const char *question, mldlg_result default_answer);
extern mldlg_result  MLRequestArgv ( MLEnvironment env, char **argv, long cardof_argv, char *buff, long size);
extern mldlg_result  MLRequestToInteract ( MLEnvironment env, mldlg_result wait_for_permission);
extern int  MLSetDialogFunction ( MLEnvironment env, long funcnum, MLDialogFunctionType func);
extern MLDialogProcPtr  MLAlertCast ( MLAlertProcPtr f);
extern MLDialogProcPtr  MLRequestCast ( MLRequestProcPtr f);
extern MLDialogProcPtr  MLConfirmCast ( MLConfirmProcPtr f);
extern MLDialogProcPtr  MLRequestArgvCast ( MLRequestArgvProcPtr f);
extern MLDialogProcPtr  MLRequestToInteractCast ( MLRequestToInteractProcPtr f);
typedef struct _mltimeval{
	unsigned long tv_sec;
	unsigned long tv_usec;
} mltimeval;
extern int  MLReady ( MLINK mlp);
typedef void *MLREADYPARALLELENV;
extern int  MLReadyParallel (MLENV, MLINK *, int, mltimeval);
typedef int (  * MLLinkWaitCallBackObject) (MLINK, void *);
extern int  MLWaitForLinkActivity (MLINK mlp);
extern int  MLWaitForLinkActivityWithCallback (MLINK mlp, MLLinkWaitCallBackObject callback);
]]
