#include "trbnet.h"
#include <thread>
#include <mutex>
#include <random>


namespace Ttrb {
	const int NOFCOMTRIES = 100;
	const int MAXFAILDELAY = 1000;
	const int MINFAILDELAY = 10;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> faildelay(MINFAILDELAY, MAXFAILDELAY);

	std::mutex TRBAccessMutex;

	inline int read_uid(uint16_t _BoardAddress, uint32_t* _buffer, size_t _buffer_size){
		int ret=-1;
		if(_buffer_size==0) {
			return ret;
		}
		for(int tries=0;
				tries<NOFCOMTRIES;
				++tries
			){
			TRBAccessMutex.lock();
			ret=trb_read_uid(_BoardAddress, _buffer, _buffer_size);
			TRBAccessMutex.unlock();
			if(ret>=0) break;
			std::cout << "FAILED COM Board Addr.: 0x" << std::hex << _BoardAddress << " trb_read_uid() " << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(faildelay(gen)));
		}
		return ret;
	}

	inline int register_write_mem(uint16_t _BoardAddress, uint16_t _RegisterAddress, uint8_t _option, uint32_t* _buffer, uint32_t _buffer_size){
		int ret=-1;
		if(_buffer_size==0) {
			return ret;
		}
		for(int tries=0;
				tries<NOFCOMTRIES;
				++tries
			){
			TRBAccessMutex.lock();
			ret=trb_register_write_mem(_BoardAddress, _RegisterAddress, _option, _buffer, _buffer_size);
			TRBAccessMutex.unlock();
			if(ret>=0) break;
			std::cout << "FAILED COM Board Addr.: 0x" << std::hex << _BoardAddress << " Register Addr.: 0x" << std::hex << _RegisterAddress << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(faildelay(gen)));
		}
		return ret;
	}

	inline int register_read_mem(uint16_t _BoardAddress, uint16_t _RegisterAddress, uint8_t _option, uint16_t _size, uint32_t* _buffer, uint32_t _buffer_size){
		int ret=-1;
		if(_buffer_size==0) {
			return ret;
		}
		for(int tries=0;
				tries<NOFCOMTRIES;
				++tries
			){
			TRBAccessMutex.lock();
			ret=trb_register_read_mem(_BoardAddress, _RegisterAddress, _option, _size, _buffer, _buffer_size);
			TRBAccessMutex.unlock();
			if(ret>=0) break;
			std::cout << "FAILED COM Board Addr.: 0x" << std::hex << _BoardAddress << " Register Addr.: 0x" << std::hex << _RegisterAddress << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(faildelay(gen)));
		}
		return ret;
	}
		
	inline int register_read(uint16_t _BoardAddress, uint16_t _RegisterAddress, uint32_t* _buffer, size_t _buffer_size){
		int ret=-1;
		if(_buffer_size==0) {
			return ret;
		}
		for(int tries=0;
				tries<NOFCOMTRIES;
				++tries
			){
			TRBAccessMutex.lock();
			ret=trb_register_read(_BoardAddress, _RegisterAddress, _buffer, _buffer_size);
			TRBAccessMutex.unlock();
			if(ret>=0) break;
			std::cout << "FAILED COM Board Addr.: 0x" << std::hex << _BoardAddress << " Register Addr.: 0x" << std::hex << _RegisterAddress << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(faildelay(gen)));
		}
		return ret;
	}

	inline int register_write(uint16_t _BoardAddress, uint16_t _RegisterAddress, uint32_t _val){
		int ret=-1;
		for(int tries=0;
				tries<NOFCOMTRIES;
				++tries
			){
			TRBAccessMutex.lock();
			ret=trb_register_write(_BoardAddress, _RegisterAddress, _val);
			TRBAccessMutex.unlock();
			if(ret!=-1) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(faildelay(gen)));
		}
		return ret;
	}
}